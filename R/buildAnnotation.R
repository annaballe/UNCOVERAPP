#' Annotate All Low Coverage Genomic Positions
#'
#' @description
#' This function identifies ALL genomic positions with low coverage (genome-wide)
#' and annotates them with variant information from a reference annotation database.
#' Unlike annotate_variants(), this processes the entire genome, not just a specific gene.
#'
#' @param sample_data Character. Path to input coverage data file in TSV format.
#' @param target_sample Character. Name of the target sample/coverage column.
#' @param coverage_threshold Numeric. Maximum coverage threshold. Default is 20.
#' @param genome Character. Genome build: "hg19" or "hg38". Default is "hg19".
#' @param annotation_file Character or NULL. Path to annotation BED file (.bed.gz).
#' @param output_intersect Character. Path for output TSV file. Default is "annotated_all_lowcov.tsv".
#' @param output_formatted Character. Path for output Excel file. Default is "annotated_all_lowcov.xlsx".
#'
#' @details
#' This function:
#' 1. Loads coverage data from input TSV file
#' 2. Filters ALL positions by coverage threshold (genome-wide)
#' 3. Queries annotation database for ALL low-coverage positions
#' 4. Annotates with functional predictions
#' 5. Generates formatted output files
#'
#' Unlike annotate_variants(), this does NOT:
#' - Query UCSC for specific gene coordinates
#' - Filter by a specific genomic region
#' - Require a gene symbol parameter
#'
#' @return Invisibly returns annotated variants data.frame
#' @export
#'
annotate_all_lowcov <- function(sample_data,
                               target_sample,
                               coverage_threshold = 20,
                               genome = "hg19",
                               annotation_file = NULL,
                               output_intersect = "annotated_all_lowcov.tsv",
                               output_formatted = "annotated_all_lowcov.xlsx") {
                                script_start <- Sys.time()

  
  # ==============================================================================
  # INPUT VALIDATION
  # ==============================================================================
  
  if (missing(sample_data)) stop("sample_data must be supplied.\n")
  if (missing(target_sample)) stop("target_sample must be supplied.\n")
  
  if (!genome %in% c("hg19", "hg38")) {
    stop("Genome must be 'hg19' or 'hg38'")
  }
  
  # Load required packages
  suppressPackageStartupMessages({
    require(dplyr)
    require(GenomicRanges)
    require(IRanges)
    require(Rsamtools)
    require(S4Vectors)
    require(openxlsx)
  })
  
  # ==============================================================================
  # CONFIGURATION
  # ==============================================================================
  
  cat("\n=== ANNOTATING ALL LOW COVERAGE POSITIONS (GENOME-WIDE) ===\n")
  cat("Sample data:", sample_data, "\n")
  cat("Target sample:", target_sample, "\n")
  cat("Coverage threshold:", coverage_threshold, "\n")
  cat("Genome:", genome, "\n")
  cat("Output intersect:", output_intersect, "\n")
  cat("Output formatted:", output_formatted, "\n\n")
  
  # ==============================================================================
  # STEP 1: LOAD INPUT DATA
  # ==============================================================================
  
  cat("=== STEP 1: LOADING INPUT DATA ===\n")
  start_time <- Sys.time()
  
  if (!file.exists(sample_data)) {
    stop("Input file not found: ", sample_data)
  }
  
  df <- tryCatch({
    read.table(sample_data, header = TRUE, sep = "\t", 
               stringsAsFactors = FALSE, check.names = FALSE)
  }, error = function(e) {
    stop("Error reading input file: ", e$message)
  })
  
  cat("Loaded", nrow(df), "rows x", ncol(df), "columns\n")
  cat("Column names:", paste(colnames(df), collapse = ", "), "\n")
  
  cat(paste("Input loading took:", 
            round(difftime(Sys.time(), start_time, units = "secs"), 2), 
            "seconds\n\n"))
  
  # ==============================================================================
  # STEP 2: FIND ACTUAL COLUMN NAME
  # ==============================================================================
  
  cat("=== STEP 2: IDENTIFYING TARGET COLUMN ===\n")
  
  actual_column <- NULL
  
  if (target_sample %in% colnames(df)) {
    actual_column <- target_sample
    cat("Found exact match for target sample:", actual_column, "\n")
  }
  
  if (is.null(actual_column)) {
    with_prefix <- paste0("count_", target_sample)
    if (with_prefix %in% colnames(df)) {
      actual_column <- with_prefix
      cat("Found target sample with 'count_' prefix:", actual_column, "\n")
    }
  }
  
  if (is.null(actual_column)) {
    without_prefix <- sub("^count_", "", target_sample)
    if (without_prefix %in% colnames(df)) {
      actual_column <- without_prefix
      cat("Found target sample after removing 'count_' prefix:", actual_column, "\n")
    }
  }
  
  if (is.null(actual_column)) {
    pattern <- gsub("^count_", "", target_sample, ignore.case = TRUE)
    matches <- grep(pattern, colnames(df), ignore.case = TRUE, value = TRUE)
    
    if (length(matches) > 0) {
      actual_column <- matches[1]
      cat("Found partial match for target sample:", actual_column, "\n")
      if (length(matches) > 1) {
        cat("WARNING: Multiple matches found:", paste(matches, collapse = ", "), "\n")
        cat("Using first match:", actual_column, "\n")
      }
    }
  }
  
  if (is.null(actual_column)) {
    stop("Target sample '", target_sample, "' not found in data.\n",
         "Available columns: ", paste(colnames(df), collapse = ", "))
  }
  
  cat("Using column:", actual_column, "\n\n")
  
  # ==============================================================================
  # STEP 3: PREPARE COVERAGE DATA (ALL POSITIONS - NO FILTERING BY REGION)
  # ==============================================================================
  
  cat("=== STEP 3: PREPARING COVERAGE DATA (ALL LOW-COVERAGE POSITIONS) ===\n")
  prep_time <- Sys.time()
  
  if (!"start" %in% colnames(df)) {
    stop("'start' column not found in input data")
  }
  if (!"chromosome" %in% colnames(df)) {
    stop("'chromosome' column not found in input data")
  }
  
  # Create bedA with all data
  bedA <- data.frame(
    chromosome = df$chromosome,
    start = df$start,
    end = df$start,
    coverage = df[[actual_column]],
    stringsAsFactors = FALSE
  )
  
  cat("Total positions in data:", nrow(bedA), "\n")
  cat("Coverage range:", min(bedA$coverage, na.rm = TRUE), "-", 
      max(bedA$coverage, na.rm = TRUE), "\n")
  cat("Chromosomes in data:", paste(unique(bedA$chromosome), collapse = ", "), "\n")
  
  bedA <- bedA %>%
    dplyr::filter(coverage <= coverage_threshold)
  
  cat("After coverage filter (<=", coverage_threshold, "):", nrow(bedA), "positions\n")
  
  if (nrow(bedA) == 0) {
    stop("No positions with coverage <= ", coverage_threshold)
  }
  
  cat("Sample of filtered data:\n")
  print(head(bedA, 10))
  
  cat(paste("\nData preparation took:",
            round(difftime(Sys.time(), prep_time, units = "secs"), 2),
            "seconds\n\n"))
  
  # ==============================================================================
  # STEP 4: SELECT ANNOTATION FILE
  # ==============================================================================
  
  cat("=== STEP 4: SELECTING ANNOTATION FILE ===\n")
  
  if (!is.null(annotation_file) && file.exists(annotation_file)) {
    file.name <- annotation_file
    cat("Using annotation file:", file.name, "\n")
  } else {
    test_hg19 <- "/home/anna/uncoverappLib/traials/geneApp/R/sorted_hg19.bed.gz"
    test_hg38 <- "/home/anna/uncoverappLib/traials/geneApp/R/sorted_hg38.bed.gz"
    
    env_hg19 <- Sys.getenv("UNCOVERAPP_HG19_ANNOTATION", unset = "")
    env_hg38 <- Sys.getenv("UNCOVERAPP_HG38_ANNOTATION", unset = "")
    
    file.name <- NULL
    
    if (identical(genome, "hg19") && nzchar(env_hg19) && file.exists(env_hg19)) {
      file.name <- env_hg19
    } else if (identical(genome, "hg38") && nzchar(env_hg38) && file.exists(env_hg38)) {
      file.name <- env_hg38
    } else if (identical(genome, "hg19") && file.exists(test_hg19)) {
      file.name <- test_hg19
    } else if (identical(genome, "hg38") && file.exists(test_hg38)) {
      file.name <- test_hg38
    }
    
    if (is.null(file.name)) {
      stop("Could not find annotation file. Please specify with annotation_file parameter")
    }
    
    cat("Using annotation file:", file.name, "\n")
  }
  
  if (!file.exists(file.name)) {
    stop("Annotation file not found: ", file.name)
  }
  
  index_file <- paste0(file.name, ".tbi")
  if (!file.exists(index_file)) {
    stop("Index file not found: ", index_file)
  }
  
  cat("Index file found:", index_file, "\n\n")
  
  # ==============================================================================
  # STEP 5: QUERY TABIX (ALL POSITIONS)
  # ==============================================================================
  
  cat("=== STEP 5: QUERYING TABIX FOR ALL LOW-COVERAGE POSITIONS ===\n")
  tabix_time <- Sys.time()
  
  bedA_for_tabix <- bedA %>%
    dplyr::mutate(chromosome = sub("^chr", "", chromosome))
  
  bedA_gr <- tryCatch({
    gr <- GenomicRanges::makeGRangesFromDataFrame(
      bedA_for_tabix,
      seqnames.field = "chromosome",
      start.field = "start",
      end.field = "end",
      keep.extra.columns = FALSE,
      ignore.strand = TRUE
    )
    GenomeInfoDb::seqlevelsStyle(gr) <- "NCBI"
    gr
  }, error = function(e) {
    stop("ERROR creating GRanges: ", e$message)
  })
  
  # Optimize query: merge nearby intervals
  cat("Optimizing: merging", length(bedA_gr), "intervals...\n")
  reduced <- GenomicRanges::reduce(bedA_gr, min.gapwidth = 100)
  cat("Reduced to", length(reduced), "intervals for query\n\n")
  
  result <- try(Rsamtools::scanTabix(file.name, param = reduced), silent = FALSE)
  
  if (inherits(result, "try-error")) {
    stop("ERROR in Tabix query")
  }
  
  lengths_result <- sapply(result, length)
  cat("Variants found:", sum(lengths_result), "\n")
  
  cat(paste("Tabix query took:",
            round(difftime(Sys.time(), tabix_time, units = "secs"), 2),
            "seconds\n\n"))
  
  if (sum(lengths_result) == 0) {
    cat("WARNING: No variants found in annotation database\n")
    write.table(data.frame(), output_intersect, 
                sep = "\t", quote = FALSE, row.names = FALSE)
    return(invisible(data.frame()))
  }
  
  # ==============================================================================
  # STEP 6: PARSE TABIX RESULTS
  # ==============================================================================
  
  cat("=== STEP 6: PARSING ANNOTATION DATA ===\n")
  parse_time <- Sys.time()
  
  dff <- lapply(result, function(elt) {
    if (length(elt) == 0) return(data.frame())
    read.csv(textConnection(elt), sep = "\t", header = FALSE, 
             stringsAsFactors = FALSE)
  })
  
  valid_dfs <- dff[sapply(dff, nrow) > 0]
  if (length(valid_dfs) == 0) {
    stop("No valid annotation data")
  }
  
  bedB <- do.call(rbind, valid_dfs)
  cat("Combined annotation:", nrow(bedB), "variants\n")
  
  # ==============================================================================
  # STEP 7: PROCESS ANNOTATION
  # ==============================================================================
  
  ncols <- ncol(bedB)
  cat("Columns in annotation:", ncols, "\n")
  
  if (ncols == 19) {
    colnames(bedB) <- c('Chromo', 'start', 'end', 'REF', 'ALT',
                        'dbsnp', 'GENENAME', 'PROTEIN_ensembl', 'field9',
                        'MutationAssessor', 'SIFT', 'Polyphen2',
                        'M_CAP', 'CADD_PHED', 'AF_gnomAD', 'ClinVar',
                        'clinvar_MedGen_id', 'HGVSc_VEP', 'HGVSp_VEP')
  } else {
    colnames(bedB) <- paste0("V", 1:ncols)
  }
  
  bedB$Chromosome <- paste0("chr", bedB[[1]])
  bedB <- bedB[, -1]
  bedB$Chromosome <- as.character(bedB$Chromosome)
  
  if ("AF_gnomAD" %in% colnames(bedB)) {
    bedB$AF_gnomAD <- suppressWarnings(as.numeric(bedB$AF_gnomAD))
  }
  if ("CADD_PHED" %in% colnames(bedB)) {
    bedB$CADD_PHED <- suppressWarnings(as.numeric(bedB$CADD_PHED))
  }
  
  cat(paste("Annotation processing took:",
            round(difftime(Sys.time(), parse_time, units = "secs"), 2),
            "seconds\n\n"))
  
  # ==============================================================================
  # STEP 8: COMPUTE OVERLAP
  # ==============================================================================
  
  cat("=== STEP 8: COMPUTING OVERLAP ===\n")
  overlap_time <- Sys.time()
  
  bed1_gr <- GenomicRanges::makeGRangesFromDataFrame(
    bedA, seqnames.field = "chromosome", ignore.strand = TRUE, 
    keep.extra.columns = TRUE
  )
  
  bed2_gr <- GenomicRanges::makeGRangesFromDataFrame(
    bedB, seqnames.field = "Chromosome", ignore.strand = TRUE, 
    keep.extra.columns = TRUE
  )
  
  tp <- GenomicRanges::findOverlaps(query = bed2_gr, subject = bed1_gr, type = "any")
  
  cat("Overlaps found:", length(tp), "\n")
  
  if (length(tp) == 0) {
    cat("WARNING: No variants overlap with low coverage regions\n")
    write.table(data.frame(), output_intersect, 
                sep = "\t", quote = FALSE, row.names = FALSE)
    return(invisible(data.frame()))
  }
  
  hits_bed2_idx <- S4Vectors::queryHits(tp)
  hits_bed1_idx <- S4Vectors::subjectHits(tp)
  
  bed1_hits <- bedA[hits_bed1_idx, ]
  bed2_hits <- bedB[hits_bed2_idx, ]
  
  intersect_df <- data.frame(
    seqnames = bed2_hits$Chromosome,
    start = bed2_hits$start,
    end = bed2_hits$end,
    coverage = bed1_hits$coverage,
    stringsAsFactors = FALSE
  )
  
  for (col in colnames(bedB)) {
    if (!col %in% c("Chromosome", "start", "end")) {
      intersect_df[[col]] <- bed2_hits[[col]]
    }
  }
  
  if (all(c("MutationAssessor", "ClinVar", "AF_gnomAD") %in% colnames(intersect_df))) {
    intersect_df$highlight_important <- grepl("H|M", intersect_df$MutationAssessor) & 
      intersect_df$ClinVar != "." & 
      !is.na(intersect_df$AF_gnomAD) & 
      intersect_df$AF_gnomAD < 0.5
    
    cat("Important variants:", sum(intersect_df$highlight_important, na.rm = TRUE), "\n")
  }
  
  cat(paste("Overlap computation took:",
            round(difftime(Sys.time(), overlap_time, units = "secs"), 2),
            "seconds\n\n"))
  
  cat("=== RESULTS ===\n")
  cat("Total annotated low-coverage positions:", nrow(intersect_df), "\n\n")
  
  # ==============================================================================
  # STEP 9: SAVE OUTPUT FILES
  # ==============================================================================
  
  cat("=== STEP 9: SAVING OUTPUT ===\n")
  
  write.table(intersect_df, output_intersect, 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  cat("Saved intersect output:", output_intersect, "\n")
  
  if (grepl("\\.xlsx?$", output_formatted, ignore.case = TRUE)) {
    cat("Creating formatted Excel file...\n")
    excel_time <- Sys.time()
    
    wb <- createWorkbook()
    addWorksheet(wb, "Low Coverage Variants")
    writeData(wb, "Low Coverage Variants", intersect_df)
    
    headerStyle <- createStyle(fontColour = "#FFFFFF", fgFill = "#4F81BD",
                               halign = "center", textDecoration = "Bold",
                               border = "TopBottomLeftRight")
    addStyle(wb, "Low Coverage Variants", headerStyle, rows = 1, 
             cols = 1:ncol(intersect_df), gridExpand = TRUE)

    setColWidths(wb, "Low Coverage Variants", cols = 1:ncol(intersect_df), widths = "auto")
    freezePane(wb, "Low Coverage Variants", firstRow = TRUE)

    # ClinVar coloring
    if ("ClinVar" %in% colnames(intersect_df)) {
      clinvar_col <- which(colnames(intersect_df) == "ClinVar")
      for (i in 1:nrow(intersect_df)) {
        style <- if (intersect_df$ClinVar[i] == ".") {
          createStyle(fgFill = "#90EE90")
        } else {
          createStyle(fgFill = "#FFB6C1")
        }
        addStyle(wb, "Low Coverage Variants", style, rows = i + 1, cols = clinvar_col)
      }
    }
    
    # CADD_PHED coloring
    if ("CADD_PHED" %in% colnames(intersect_df)) {
      cadd_col <- which(colnames(intersect_df) == "CADD_PHED")
      for (i in 1:nrow(intersect_df)) {
        val <- intersect_df$CADD_PHED[i]
        style <- if (!is.na(val) && val > 20) {
          createStyle(fgFill = "#FFB6C1")
        } else if (!is.na(val)) {
          createStyle(fgFill = "#90EE90")
        } else {
          NULL
        }
        if (!is.null(style)) addStyle(wb, "Low Coverage Variants", style, rows = i + 1, cols = cadd_col)
      }
    }
    
    saveWorkbook(wb, output_formatted, overwrite = TRUE)
    cat("Saved formatted output:", output_formatted, "\n")
    cat(paste("Excel formatting took:", 
              round(difftime(Sys.time(), excel_time, units = "secs"), 2), 
              "seconds\n"))
  }
  cat("Total script execution:",round(difftime(Sys.time(), script_start, units = "secs"), 2),"seconds\n")
  cat("\n=== DONE ===\n")
  invisible(intersect_df)
}
