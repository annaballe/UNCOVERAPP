filtered_low_nucl <- reactive({
  if (is.null(mysample()))
    return(NULL)
  
  df <- mysample()
  thr <- input$coverage_co
  chr_val <- input$Chromosome
  
  cat("Filtering for chromosome:", chr_val, "with coverage threshold:", thr, "\n")
  cat("Sample data dimensions:", dim(df), "\n")
  cat("Available chromosomes in data:", unique(df$chromosome), "\n")
  
  if (identical(thr, "all")) {
    result <- dplyr::filter(df, chromosome == chr_val)
  } else {
    result <- dplyr::filter(df, chromosome == chr_val,
                           coverage <= as.numeric(thr))
  }
  
  cat("Filtered result dimensions:", dim(result), "\n")
  if (nrow(result) > 0) {
    cat("Sample of filtered data:\n")
    print(head(result))
  } else {
    cat("No data found for chromosome", chr_val, "with coverage <=", thr, "\n")
  }
  
  return(result)
})


# ============================================================================
# REACTIVE BASE: Calcolo una sola volta
# ============================================================================

# ============================================================================
# compute-annotation.R - ANNOTAZIONE VARIANTI
# ============================================================================

# ============================================================================
# STEP 1: QUERY ANNOTATION DATABASE (calcola UNA SOLA VOLTA)
# ============================================================================

intBED <- reactive({
  cat("\n=== ANNOTATION: intBED (ONE-TIME CALCULATION) ===\n")
  
  if (is.null(filtered_low_nucl())) {
    cat("ERROR: No filtered data\n")
    return(NULL)
  }
  
  bedA <- filtered_low_nucl()
  cat("Low coverage positions:", nrow(bedA), "\n")
  
  if (nrow(bedA) == 0) {
    cat("ERROR: No positions to annotate\n")
    return(NULL)
  }
  
  # FILTER BY USER REGION (optional)
  query <- input$query_Database
  
  if (!is.null(query) && query != "") {
    cat("Filtering by region:", query, "\n")
    
    query.regions <- tryCatch({
      read.table(text=gsub("[:-]+", " ", query, perl=TRUE),
                 header=FALSE, col.names = c("chr", "start", "end"))
    }, error = function(e) NULL)
    
    if (!is.null(query.regions)) {
      query_chr <- if (grepl("^chr", query.regions$chr)) {
        query.regions$chr
      } else {
        paste0("chr", query.regions$chr)
      }
      
      bedA <- bedA %>%
        dplyr::filter(
          chromosome == query_chr,
          end >= query.regions$start,
          start <= query.regions$end
        )
      
      if (nrow(bedA) == 0) {
        cat("WARNING: No data in region\n")
        return(NULL)
      }
    }
  }
  
  # SELECT ANNOTATION FILE
  test_hg19 <- "/home/anna/uncoverappLib/traials/geneApp/R/sorted_hg19.bed.gz"
  test_hg38 <- "/home/anna/uncoverappLib/traials/geneApp/R/sorted_hg38.bed.gz"
  env_hg19 <- Sys.getenv("UNCOVERAPP_HG19_ANNOTATION", unset = "")
  env_hg38 <- Sys.getenv("UNCOVERAPP_HG38_ANNOTATION", unset = "")
  
  file.name <- NULL
  
  if (identical(input$UCSC_Genome, "hg19") && nzchar(env_hg19) && file.exists(env_hg19)) {
    file.name <- env_hg19
  } else if (identical(input$UCSC_Genome, "hg38") && nzchar(env_hg38) && file.exists(env_hg38)) {
    file.name <- env_hg38
  } else if (identical(input$UCSC_Genome, "hg19") && file.exists(test_hg19)) {
    file.name <- test_hg19
  } else if (identical(input$UCSC_Genome, "hg38") && file.exists(test_hg38)) {
    file.name <- test_hg38
  } else {
    m <- uncoverappLib::getAnnotationFiles()
    data_files <- m[grepl("\\.gz$", m) & !grepl("\\.tbi$", m)]
    
    if (length(data_files) > 0) {
      idx <- if (input$UCSC_Genome == "hg19") {
        grep("hg19|GRCh37", data_files, ignore.case = TRUE)
      } else {
        grep("hg38|GRCh38", data_files, ignore.case = TRUE)
      }
      file.name <- data_files[if(length(idx) > 0) idx[1] else 1]
    }
  }
  
  cat("Annotation file:", file.name, "\n")
  
  # QUERY TABIX
  cat("Querying Tabix...\n")
  
  bedA_for_tabix <- bedA %>%
    dplyr::mutate(chromosome = sub("^chr", "", chromosome))
  
  bedA_gr <- tryCatch({
    GenomicRanges::makeGRangesFromDataFrame(
      bedA_for_tabix,
      seqnames.field = "chromosome",
      start.field = "start",
      end.field = "end",
      keep.extra.columns = FALSE,
      ignore.strand = TRUE
    )
  }, error = function(e) {
    cat("ERROR creating GRanges:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(bedA_gr)) return(NULL)
  
  query_gr <- if (length(bedA_gr) > 1000) {
    cat("Optimizing: merging intervals...\n")
    reduced <- GenomicRanges::reduce(bedA_gr, min.gapwidth = 100)
    cat("Reduced from", length(bedA_gr), "to", length(reduced), "\n")
    reduced
  } else {
    bedA_gr
  }
  
  result <- try(Rsamtools::scanTabix(file.name, param = query_gr), silent = FALSE)
  
  if (inherits(result, "try-error")) {
    cat("ERROR in Tabix\n")
    return(NULL)
  }
  
  lengths_result <- sapply(result, length)
  cat("Variants found:", sum(lengths_result), "\n")
  
  if (sum(lengths_result) == 0) {
    cat("No variants\n")
    return(NULL)
  }
  
  dff <- lapply(result, function(elt) {
    if (length(elt) == 0) return(data.frame())
    read.csv(textConnection(elt), sep="\t", header=FALSE, stringsAsFactors = FALSE)
  })
  
  valid_dfs <- dff[sapply(dff, nrow) > 0]
  if (length(valid_dfs) == 0) return(NULL)
  
  bedB <- do.call(rbind, valid_dfs)
  cat("Combined annotation:", nrow(bedB), "variants\n")
  
  # PROCESS ANNOTATION
  ncols <- ncol(bedB)
  
  if (ncols == 19) {
    colnames(bedB) <- c('Chromo', 'start','end','REF','ALT',
                        'dbsnp','GENENAME', 'PROTEIN_ensembl',
                        'MutationAssessor','SIFT','Polyphen2',
                        'M_CAP','CADD_PHED','AF_gnomAD','ClinVar',
                        'clinvar_MedGen_id','clinvar_OMIM_id','HGVSc_VEP','HGVSp_VEP')
  } else {
    colnames(bedB) <- paste0("V", 1:ncols)
  }
  
  bedB$Chromosome <- paste0("chr", bedB[[1]])
  bedB <- bedB[, -1]
  bedB$Chromosome <- as.character(bedB$Chromosome)
  bedB$AF_gnomAD <- suppressWarnings(as.numeric(bedB$AF_gnomAD))
  bedB$CADD_PHED <- suppressWarnings(as.numeric(bedB$CADD_PHED))
  
  # OVERLAP
  cat("Computing overlap...\n")
  
  bed1_gr <- GenomicRanges::makeGRangesFromDataFrame(
    bedA, seqnames.field = "chromosome", ignore.strand = TRUE, keep.extra.columns = TRUE
  )
  
  bed2_gr <- GenomicRanges::makeGRangesFromDataFrame(
    bedB, seqnames.field = "Chromosome", ignore.strand = TRUE, keep.extra.columns = TRUE
  )
  
  tp <- GenomicRanges::findOverlaps(query = bed2_gr, subject = bed1_gr, type = "any")
  
  cat("Overlaps:", length(tp), "\n")
  
  if (length(tp) == 0) return(data.frame())
  
  hits_bed2_idx <- S4Vectors::queryHits(tp)
  hits_bed1_idx <- S4Vectors::subjectHits(tp)
  
  bed1_hits <- bedA[hits_bed1_idx, ]
  bed2_hits <- bedB[hits_bed2_idx, ]
  
  intersect_df <- data.frame(
    seqnames = bed2_hits$Chromosome,
    start = bed2_hits$start,
    end = bed2_hits$end,
    coverage = bed1_hits$coverage,
    REF = bed2_hits$REF,
    ALT = bed2_hits$ALT,
    dbsnp = bed2_hits$dbsnp,
    GENENAME = bed2_hits$GENENAME,
    PROTEIN_ensembl = bed2_hits$PROTEIN_ensembl,
    MutationAssessor = bed2_hits$MutationAssessor,
    SIFT = bed2_hits$SIFT,
    Polyphen2 = bed2_hits$Polyphen2,
    M_CAP = bed2_hits$M_CAP,
    CADD_PHED = bed2_hits$CADD_PHED,
    AF_gnomAD = bed2_hits$AF_gnomAD,
    ClinVar = bed2_hits$ClinVar,
    clinvar_MedGen_id = bed2_hits$clinvar_MedGen_id,
    HGVSc_VEP = bed2_hits$HGVSc_VEP,
    HGVSp_VEP = bed2_hits$HGVSp_VEP,
    stringsAsFactors = FALSE
  )
  
  cat("Final result:", nrow(intersect_df), "variants\n")
  return(intersect_df)
})

# ============================================================================
# STEP 2: CLEAN DATAFRAME (per operazioni dplyr, export, ecc.)
# ============================================================================

annotated_variants_data <- reactive({
  cat("\n=== ANNOTATION: Preparing clean dataframe ===\n")
  
  data <- intBED()
  
  if (is.null(data) || nrow(data) == 0) {
    return(NULL)
  }
  
  # Aggiungi colonna helper per highlight
  data$highlight_important <- grepl("H|M", data$MutationAssessor) & 
                              data$ClinVar != "." & 
                              !is.na(data$AF_gnomAD) & 
                              data$AF_gnomAD < 0.5
  
  cat("Dataframe ready:", nrow(data), "variants\n")
  return(data)
})

# ============================================================================
# STEP 3: DT WIDGET (SOLO per visualizzazione)
# ============================================================================

condform_table <- reactive({
  cat("\n=== ANNOTATION: Building DT widget ===\n")
  
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Formatting table...", value = 0.5)
  
  data <- annotated_variants_data()
  
  if (is.null(data) || nrow(data) == 0) {
    return(NULL)
  }
  
  cat("Formatting", nrow(data), "rows\n")
  
  dt <- DT::datatable(
    data,
    options = list(
      pageLength = 25,
      scrollX = TRUE,
      scrollY = "600px",
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel'),
      columnDefs = list(
        list(targets = which(colnames(data) == "highlight_important") - 1, visible = FALSE)
      )
    ),
    extensions = 'Buttons',
    rownames = FALSE
  ) %>%
    DT::formatStyle('ClinVar',
                    backgroundColor = DT::styleEqual(c('.'), c('lightgreen'), default = 'lightcoral')) %>%
    DT::formatStyle('CADD_PHED',
                    backgroundColor = DT::styleInterval(c(20), c('lightgreen', 'lightcoral'))) %>%
    DT::formatStyle('MutationAssessor',
                    backgroundColor = DT::styleEqual(c('H', 'M'), c('lightcoral', 'lightyellow'), default = 'lightgreen')) %>%
    DT::formatStyle('M_CAP',
                    backgroundColor = DT::styleEqual(c('D'), c('lightcoral'), default = 'lightgreen')) %>%
    DT::formatStyle('AF_gnomAD',
                    backgroundColor = DT::styleInterval(c(0.5), c('lightcoral', 'lightgreen'))) %>%
    DT::formatStyle(c('start', 'end'), 'highlight_important',
                    backgroundColor = DT::styleEqual(c(TRUE, FALSE), c('yellow', 'white'))) %>%
    DT::formatStyle(c('start', 'end'), 'highlight_important',
                    color = DT::styleEqual(c(TRUE, FALSE), c('red', 'green')))
  
  cat("DT widget created\n")
  return(dt)
})