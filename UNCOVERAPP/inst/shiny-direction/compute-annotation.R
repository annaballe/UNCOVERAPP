filtered_low_nucl<- reactive ({
  if (is.null(mysample()))
    return(NULL)
  
  df <- mysample()
  thr <- input$coverage_co
  
  # Get chromosome value - use input directly to avoid reactive timing issues
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

intBED<- reactive({
  cat("=== STARTING intBED FUNCTION ===\n")
  
  if (is.null(filtered_low_nucl())) {
    cat("ERROR: filtered_low_nucl() returned NULL\n")
    return(NULL)
  }
  
  bedA<- filtered_low_nucl()
  cat("bedA (coverage data) dimensions:", dim(bedA), "\n")
  cat("bedA sample data:\n")
  print(head(bedA))

  # TEST OVERRIDE: prefer custom local files or env vars for annotation during testing
  test_hg19 <- "/home/anna/uncoverappLib/traials/geneApp/R/sorted_hg19.bed.gz"
  test_hg38 <- "/home/anna/uncoverappLib/traials/geneApp/R/sorted_hg38.bed.gz"
  env_hg19 <- Sys.getenv("UNCOVERAPP_HG19_ANNOTATION", unset = "")
  env_hg38 <- Sys.getenv("UNCOVERAPP_HG38_ANNOTATION", unset = "")

  file.name <- NULL
  # 1) Environment variable overrides (if set)
  if (identical(input$UCSC_Genome, "hg19") && nzchar(env_hg19) && file.exists(env_hg19)) {
    file.name <- env_hg19
  }
  if (identical(input$UCSC_Genome, "hg38") && nzchar(env_hg38) && file.exists(env_hg38)) {
    file.name <- env_hg38
  }
  # 2) Test files in repo (if present)
  if (is.null(file.name) && identical(input$UCSC_Genome, "hg19") && file.exists(test_hg19)) {
    file.name <- test_hg19
  }
  if (is.null(file.name) && identical(input$UCSC_Genome, "hg38") && file.exists(test_hg38)) {
    file.name <- test_hg38
  }
  # 3) Fallback to default cache selection
  if (is.null(file.name)) {
    m <- uncoverappLib::getAnnotationFiles()
    # Separate data files (.gz) from index files (.tbi)
    data_files <- m[grepl("\\.gz$", m) & !grepl("\\.tbi$", m)]
    # Selection rules:
    # - If multiple data files exist (e.g., hg19 + hg38), pick by genome string.
    # - Else, use the only data file present (ignore genome choice).
    if (length(data_files) > 1) {
      idx <- if (input$UCSC_Genome == "hg19") {
        grep("hg19|GRCh37", data_files, ignore.case = TRUE)
      } else {
        grep("hg38|GRCh38", data_files, ignore.case = TRUE)
      }
      if (length(idx) == 0) idx <- 1
      file.name <- data_files[idx[1]]
    } else if (length(data_files) == 1) {
      file.name <- data_files[1]
    } else {
      # Fallback: if nothing matched (unexpected), try the first returned path
      file.name <- m[1]
    }
  }
  # Print to console which file is used (useful during testing)
  cat("Using annotation file:", file.name, "\n")

  #second and tirth columns are hg19 positions
  cat("\n=== PROCESSING QUERY ===\n")
  query <- c(input$query_Database)
  cat("Original query string:", query, "\n")
  
  query.regions= read.table(text=gsub("[:-]+", " ", query, perl=TRUE),
             header=FALSE, col.names = c("chr", "start", "end"))
  
  if (is.null(query.regions)) {
    cat("ERROR: Failed to parse query regions\n")
    return(NULL)
  }
  
  cat("Parsed query regions:\n")
  print(query.regions)
  cat("\n=== RUNNING TABIX QUERY ===\n")
  result<- try({
    cat("Creating GRanges from query regions...\n")
    fq= GenomicRanges::makeGRangesFromDataFrame(query.regions, keep.extra.columns = TRUE)
    cat("GRanges created successfully\n")
    
    cat("Running scanTabix with file:", file.name, "\n")
    res <- Rsamtools::scanTabix(file.name, param=fq)
    cat("Tabix completed\n")
    
    lengths_result <- sapply(res, length)
    cat("Tabix result lengths:", lengths_result, "\n")
    dff <- Map(function(elt) {
      if (length(elt) == 0) {
        return(data.frame())
      }
      read.csv(textConnection(elt), sep="\t", header=FALSE, stringsAsFactors = FALSE)
    }, res)
    
    # Check if we have any data
    valid_dfs <- dff[sapply(dff, nrow) > 0]
    if (length(valid_dfs) == 0) {
      cat("No annotation data found for query region\n")
      return(NULL)
    }
    
    bedB <- do.call(rbind, valid_dfs)

  })
  if ("try-error" %in% class(result)) {
    err_msg <- 'no coordinates recognized'
    cat("Error in annotation processing:", err_msg, "\n")
    return(NULL)
  }
  
  # Check if bedB exists and has data
  if (is.null(bedB) || nrow(bedB) == 0) {
    cat("No annotation data to process\n")
    return(NULL)
  }

  print(head(bedB))
  # Check actual number of columns and set names accordingly
  ncols <- ncol(bedB)
  cat("Number of columns in annotation data:", ncols, "\n")
  
  if (ncols == 19) {
    colnames(bedB)<- c ('Chromo', 'start','end','REF','ALT',
                        'dbsnp','GENENAME', 'PROTEIN_ensembl', 'field9',
                        'MutationAssessor','SIFT','Polyphen2',
                        'M_CAP','CADD_PHED','AF_gnomAD','ClinVar',
                        'clinvar_MedGen_id','HGVSc_VEP','HGVSp_VEP')
  } else {
    # Fallback: use generic column names
    colnames(bedB) <- paste0("V", 1:ncols)
  }
  str(bedB)
  
  # For this data structure, start and end are already in the right columns
  # No need to rename - they're already named 'start' and 'end'
  cat("Column names after assignment:", colnames(bedB), "\n")

  for (i in bedB[1]){
    Chromosome<-paste("chr", i, sep="")
    bedB<- cbind(Chromosome, bedB)
    bedB[,2]<-NULL
  }
  bedB$Chromosome= as.character(bedB$Chromosome)
  bedB$AF_gnomAD= suppressWarnings(as.numeric(bedB$AF_gnomAD))
  bedB$CADD_PHED= suppressWarnings(as.numeric(bedB$CADD_PHED))
  intersectBedFiles.GR <- function(bed1,bed2) {
    require(GenomicRanges)
    require(IRanges)
    bed1 <- makeGRangesFromDataFrame(bedA,ignore.strand = TRUE,
                                     keep.extra.columns = TRUE)
    bed2 <- makeGRangesFromDataFrame(bedB,ignore.strand = TRUE,
                                     keep.extra.columns = TRUE)
    tp<- findOverlaps(query = bed1, subject = bed2, type="any")
    intersect_df = data.frame(bed1[queryHits(tp),], bed2[subjectHits(tp),])
    return(intersect_df)
  }
  intersect_df<- intersectBedFiles.GR(bedA, bedB)
  return(intersect_df)
})

#make reactive dataframe

condform_table<- reactive ({
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "table construction in progress",
               detail = 'This may take a while', value = 0)
  Sys.sleep(0.1)
  library(condformat)
  if (is.null(intBED()))
    return(NULL)
  condformat(intBED()) %>%
    rule_fill_discrete(ClinVar, expression= ClinVar !=".",
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(CADD_PHED, expression= CADD_PHED > 20,
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(MutationAssessor, expression=  MutationAssessor =='H',
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(M_CAP, expression=  M_CAP =='D',
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(AF_gnomAD, expression=
                         ifelse(is.na(AF_gnomAD) | AF_gnomAD < 0.5,
                                'TRUE','FALSE') ,
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(c(start, end),
                       expression = grepl("H|M", MutationAssessor) &
                         ClinVar !="." &  AF_gnomAD < 0.5 ,
                       colours= c("TRUE"= "yellow", "FALSE"= ""))%>%
    rule_css(c(start, end),
             expression = ifelse(grepl("H|M", MutationAssessor) &
                                   ClinVar !="." &  AF_gnomAD < 0.5,
                                 "red", "green"),
             css_field = "color")
})
