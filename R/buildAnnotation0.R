#' Annotate Genomic Variants with Coverage and Functional Information
#'
#' @description
#' This function identifies genomic positions with low coverage and annotates them
#' with variant information from a reference annotation database. It performs
#' intersection between coverage data and variant annotations, highlighting
#' potentially important variants based on multiple criteria.
#'
#' @return Invisibly returns a data.frame with annotated variants. Also writes
#'   two output files: a TSV file with all results and a formatted Excel file.
#'
#' @param sample_data Character. Path to input coverage data file in TSV format.
#'   Must contain columns: 'chromosome', 'start', and coverage columns.
#' @param target_sample Character. Name of the target sample/chromosome column
#'   (e.g., "chr10_positions2"). The function will automatically handle variations
#'   with/without "count_" prefix.
#' @param gene Character. Gene symbol to query (e.g., "BRCA1"). Gene coordinates
#'   will be retrieved from UCSC/Ensembl databases.
#' @param coverage_threshold Numeric. Maximum coverage threshold for filtering
#'   positions. Default is 20. Positions with coverage <= threshold are retained.
#' @param query_region Character. Optional genomic region in format
#'   "chr:start-end" (e.g., "chr10:73697750-73697765"). If empty, gene
#'   coordinates will be used. Default is "".
#' @param genome Character. Genome build: "hg19" or "hg38". Default is "hg19".
#' @param annotation_file Character or NULL. Path to annotation BED file
#'   (.bed.gz with .tbi index). If NULL, function will search for default
#'   annotation files. Default is NULL.
#' @param output_intersect Character. Path for output TSV file with intersected
#'   results. Default is "intersect_output.tsv".
#' @param output_formatted Character. Path for output Excel file with formatting.
#'   Default is "formatted_output.xlsx".
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Loads coverage data from input TSV file
#'   \item Retrieves gene coordinates from UCSC/Ensembl (requires biomaRt or rtracklayer)
#'   \item Filters positions by coverage threshold and genomic region
#'   \item Queries annotation database using Tabix
#'   \item Identifies overlapping variants
#'   \item Annotates variants with functional predictions
#'   \item Highlights important variants based on multiple criteria
#'   \item Generates formatted output files
#' }
#'
#' @section Important Variants:
#' Variants are flagged as "important" if they meet ALL of these criteria:
#' \itemize{
#'   \item MutationAssessor score is "H" (high) or "M" (medium)
#'   \item ClinVar annotation is present (not ".")
#'   \item gnomAD allele frequency < 0.5
#' }
#'
#' @section Excel Formatting:
#' The formatted Excel output includes color-coded cells:
#' \itemize{
#'   \item ClinVar: Red if present, Green if absent
#'   \item CADD_PHED: Red if > 20, Green otherwise
#'   \item MutationAssessor: Red for "H", Yellow for "M", Green otherwise
#'   \item M_CAP: Red for "D" (deleterious), Green otherwise
#'   \item AF_gnomAD: Red if < 0.5 (rare), Green otherwise
#'   \item Start/End positions: Yellow for important variants
#' }
#'
#' @section Required Packages:
#' \itemize{
#'   \item dplyr
#'   \item GenomicRanges
#'   \item IRanges
#'   \item Rsamtools
#'   \item S4Vectors
#'   \item openxlsx
#'   \item biomaRt or rtracklayer (for gene coordinate retrieval)
#' }
#'
#' @section Annotation File:
#' The annotation file must be:
#' \itemize{
#'   \item Tab-separated BED format
#'   \item Bgzip compressed (.bed.gz)
#'   \item Tabix indexed (.bed.gz.tbi)
#' }
#'
#' Expected columns (19 total):
#' Chromosome, start, end, REF, ALT, dbsnp, GENENAME, PROTEIN_ensembl,
#' field9, MutationAssessor, SIFT, Polyphen2, M_CAP, CADD_PHED,
#' AF_gnomAD, ClinVar, clinvar_MedGen_id, HGVSc_VEP, HGVSp_VEP
#'
#' @export
#' @import Rsamtools
#' @import GenomicRanges
#' @import IRanges
#' @import S4Vectors
#' @import openxlsx
#' @importFrom dplyr filter mutate
#' @importFrom utils read.table write.table
#'
#' @examples
#' \dontrun{
#' # Basic usage with gene coordinates
#' annotate_variants(
#'   sample_data = "coverage.tsv",
#'   target_sample = "chr10_positions2",
#'   gene = "BRCA1",
#'   coverage_threshold = 20,
#'   genome = "hg19"
#' )
#'
#' # With custom region
#' annotate_variants(
#'   sample_data = "coverage.tsv",
#'   target_sample = "chr10_positions2",
#'   gene = "BRCA1",
#'   coverage_threshold = 15,
#'   query_region = "chr10:73697750-73697765",
#'   genome = "hg19",
#'   annotation_file = "/path/to/sorted_hg19.bed.gz",
#'   output_intersect = "results.tsv",
#'   output_formatted = "results.xlsx"
#' )
#' }
#'
#' @author Anna
#' @seealso \code{\link[GenomicRanges]{findOverlaps}}, \code{\link[Rsamtools]{scanTabix}}
annotate_variants <- function(sample_data,
                             target_sample,
                             gene,
                             coverage_threshold = 20,
                             query_region = "",
                             genome = "hg19",
                             annotation_file = NULL,
                             output_intersect = "intersect_output.tsv",
                             output_formatted = "formatted_output.xlsx") {
  
  # ==============================================================================
  # INPUT VALIDATION
  # ==============================================================================
  
  if (missing(sample_data)) stop("sample_data must be supplied.\n")
  if (missing(target_sample)) stop("target_sample must be supplied.\n")
  if (missing(gene)) stop("gene must be supplied.\n")
  
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
  
  cat("\n=== CONFIGURATION ===\n")
  cat("Sample data:", sample_data, "\n")
  cat("Target sample:", target_sample, "\n")
  cat("Gene:", gene, "\n")
  cat("Coverage threshold:", coverage_threshold, "\n")
  cat("Query region:", ifelse(query_region == "", "(none - will use gene coordinates)", query_region), "\n")
  cat("Genome:", genome, "\n")
  cat("Annotation file:", ifelse(is.null(annotation_file), "(auto-detect)", annotation_file), "\n")
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
  
  # Read coverage data
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
  # STEP 2: FIND ACTUAL COLUMN NAME (handle count_ prefix)
  # ==============================================================================
  
  cat("=== STEP 2: IDENTIFYING TARGET COLUMN ===\n")
  
  actual_column <- NULL
  
  # Option 1: Exact match
  if (target_sample %in% colnames(df)) {
    actual_column <- target_sample
    cat("Found exact match for target sample:", actual_column, "\n")
  }
  
  # Option 2: Try with count_ prefix
  if (is.null(actual_column)) {
    with_prefix <- paste0("count_", target_sample)
    if (with_prefix %in% colnames(df)) {
      actual_column <- with_prefix
      cat("Found target sample with 'count_' prefix:", actual_column, "\n")
    }
  }
  
  # Option 3: Try removing count_ prefix from user input
  if (is.null(actual_column)) {
    without_prefix <- sub("^count_", "", target_sample)
    if (without_prefix %in% colnames(df)) {
      actual_column <- without_prefix
      cat("Found target sample after removing 'count_' prefix:", actual_column, "\n")
    }
  }
  
  # Option 4: Case-insensitive partial match
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
  
  # If still not found, show error
  if (is.null(actual_column)) {
    stop("Target sample '", target_sample, "' not found in data.\n",
         "Tried variations: '", target_sample, "', 'count_", target_sample, "'\n",
         "Available columns: ", paste(colnames(df), collapse = ", "))
  }
  
  cat("Using column:", actual_column, "\n\n")
  
  # ==============================================================================
  # STEP 3: GET GENE COORDINATES FROM UCSC
  # ==============================================================================
  
  cat("=== STEP 3: QUERYING UCSC FOR GENE COORDINATES ===\n")
  coord_time <- Sys.time()
  cat("Gene:", gene, "\n")
  cat("Genome:", genome, "\n")
  
  gene_coords <- NULL
  
  # Try biomaRt first
  if (requireNamespace("biomaRt", quietly = TRUE)) {
    tryCatch({
      if (genome == "hg19") {
        mart <- biomaRt::useEnsembl(biomart = "ensembl", 
                                    dataset = "hsapiens_gene_ensembl",
                                    version = "GRCh37")
      } else {
        mart <- biomaRt::useEnsembl(biomart = "ensembl", 
                                    dataset = "hsapiens_gene_ensembl")
      }
      
      gene_info <- biomaRt::getBM(
        attributes = c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
        filters = "hgnc_symbol",
        values = gene,
        mart = mart
      )
      
      if (nrow(gene_info) > 0) {
        chr <- paste0("chr", gene_info$chromosome_name[1])
        start <- gene_info$start_position[1]
        end <- gene_info$end_position[1]
        
        cat("Found gene coordinates:\n")
        cat("  Chromosome:", chr, "\n")
        cat("  Start:", start, "\n")
        cat("  End:", end, "\n\n")
        
        gene_coords <- list(chr = chr, start = start, end = end)
      }
    }, error = function(e) {
      cat("biomaRt query failed:", e$message, "\n")
    })
  }
  
  # Fallback: try rtracklayer
  if (is.null(gene_coords) && requireNamespace("rtracklayer", quietly = TRUE)) {
    tryCatch({
      session <- rtracklayer::browserSession("UCSC")
      rtracklayer::genome(session) <- genome
      
      query <- rtracklayer::ucscTableQuery(session, "refGene")
      genes <- rtracklayer::getTable(query)
      
      gene_data <- genes[grepl(gene, genes$name2, ignore.case = TRUE), ]
      
      if (nrow(gene_data) > 0) {
        chr <- gene_data$chrom[1]
        start <- min(gene_data$txStart)
        end <- max(gene_data$txEnd)
        
        cat("Found gene coordinates:\n")
        cat("  Chromosome:", chr, "\n")
        cat("  Start:", start, "\n")
        cat("  End:", end, "\n\n")
        
        gene_coords <- list(chr = chr, start = start, end = end)
      }
    }, error = function(e) {
      cat("rtracklayer query failed:", e$message, "\n")
    })
  }
  
  if (is.null(gene_coords)) {
    warning("Could not query UCSC automatically. Please install biomaRt or rtracklayer package.")
    stop("Could not find coordinates for gene: ", gene)
  }
  
  cat(paste("Gene coordinate retrieval took:",
            round(difftime(Sys.time(), coord_time, units = "secs"), 2),
            "seconds\n\n"))
  
  # ==============================================================================
  # STEP 4: DETERMINE QUERY REGION
  # ==============================================================================
  
  if (query_region != "") {
    cat("=== USING USER-PROVIDED QUERY REGION ===\n")
    query_string <- query_region
  } else {
    cat("=== USING GENE COORDINATES AS QUERY REGION ===\n")
    query_string <- paste0(gene_coords$chr, ":", gene_coords$start, "-", gene_coords$end)
  }
  
  cat("Query region:", query_string, "\n\n")
  
  # Parse query region
  query.regions <- tryCatch({
    read.table(text = gsub("[:-]+", " ", query_string, perl = TRUE),
               header = FALSE, col.names = c("chr", "start", "end"))
  }, error = function(e) {
    stop("Failed to parse query region: ", query_string)
  })
  
  cat("Parsed query:\n")
  cat("  Chromosome:", query.regions$chr, "\n")
  cat("  Start:", query.regions$start, "\n")
  cat("  End:", query.regions$end, "\n\n")
  
  # ==============================================================================
  # STEP 5: PREPARE COVERAGE DATA
  # ==============================================================================
  
  cat("=== STEP 5: PREPARING COVERAGE DATA ===\n")
  prep_time <- Sys.time()
  
  # Validate required columns
  if (!"start" %in% colnames(df)) {
    stop("'start' column not found in input data")
  }
  if (!"chromosome" %in% colnames(df)) {
    stop("'chromosome' column not found in input data")
  }
  
  # Ensure chromosome format matches
  query_chr <- if (grepl("^chr", query.regions$chr)) {
    query.regions$chr
  } else {
    paste0("chr", query.regions$chr)
  }
  
  # Create bedA: chromosome from DATA, start, end, coverage
  bedA <- data.frame(
    chromosome = df$chromosome,
    start = df$start,
    end = df$start,  # For single positions, start == end
    coverage = df[[actual_column]],
    stringsAsFactors = FALSE
  )
  
  cat("Coverage data dimensions:", nrow(bedA), "rows\n")
  cat("Coverage range:", min(bedA$coverage, na.rm = TRUE), "-", 
      max(bedA$coverage, na.rm = TRUE), "\n")
  cat("Chromosomes in data:", paste(unique(bedA$chromosome), collapse = ", "), "\n")
  cat("Query chromosome:", query_chr, "\n")
  
  # Filter by coverage threshold
  bedA <- bedA %>%
    dplyr::filter(coverage <= coverage_threshold)
  
  cat("After coverage filter (<=", coverage_threshold, "):", nrow(bedA), "positions\n")
  
  if (nrow(bedA) == 0) {
    stop("No positions with coverage <= ", coverage_threshold)
  }
  
  # Filter by query region
  bedA <- bedA %>%
    dplyr::filter(
      chromosome == query_chr,
      end >= query.regions$start,
      start <= query.regions$end
    )
  
  cat("After region filter:", nrow(bedA), "positions\n")
  
  if (nrow(bedA) == 0) {
    stop("No positions in query region with coverage <= ", coverage_threshold)
  }
  
  cat("Sample of filtered data:\n")
  print(head(bedA))
  
  cat(paste("\nData preparation took:",
            round(difftime(Sys.time(), prep_time, units = "secs"), 2),
            "seconds\n\n"))
  
  # ==============================================================================
  # STEP 6: SELECT ANNOTATION FILE
  # ==============================================================================
  
  cat("=== STEP 6: SELECTING ANNOTATION FILE ===\n")
  
  if (!is.null(annotation_file) && file.exists(annotation_file)) {
    file.name <- annotation_file
    cat("Using annotation file:", file.name, "\n")
  } else {
    # Try default locations
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
    } else if (requireNamespace("uncoverappLib", quietly = TRUE)) {
      m <- uncoverappLib::getAnnotationFiles()
      data_files <- m[grepl("\\.gz$", m) & !grepl("\\.tbi$", m)]
      
      if (length(data_files) > 0) {
        idx <- if (genome == "hg19") {
          grep("hg19|GRCh37", data_files, ignore.case = TRUE)
        } else {
          grep("hg38|GRCh38", data_files, ignore.case = TRUE)
        }
        file.name <- data_files[if(length(idx) > 0) idx[1] else 1]
      }
    }
    
    if (is.null(file.name)) {
      stop("Could not find annotation file. Please specify with annotation_file parameter")
    }
    
    cat("Using annotation file:", file.name, "\n")
  }
  
  if (!file.exists(file.name)) {
    stop("Annotation file not found: ", file.name)
  }
  
  # Check for index file
  index_file <- paste0(file.name, ".tbi")
  if (!file.exists(index_file)) {
    stop("Index file not found: ", index_file, 
         "\nPlease create with: tabix -p bed ", file.name)
  }
  
  cat("Index file found:", index_file, "\n\n")
  
  # ==============================================================================
  # STEP 7: QUERY TABIX
  # ==============================================================================
  
  cat("=== STEP 7: QUERYING TABIX ===\n")
  tabix_time <- Sys.time()
  
  # Remove "chr" prefix for tabix query
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
  
  # Optimize query
  query_gr <- if (length(bedA_gr) > 1000) {
    cat("Optimizing: merging", length(bedA_gr), "intervals...\n")
    reduced <- GenomicRanges::reduce(bedA_gr, min.gapwidth = 100)
    cat("Reduced to", length(reduced), "intervals\n")
    reduced
  } else {
    bedA_gr
  }
  
  result <- try(Rsamtools::scanTabix(file.name, param = query_gr), silent = FALSE)
  
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
    cat("Created empty output file:", output_intersect, "\n")
    return(invisible(data.frame()))
  }
  
  # ==============================================================================
  # STEP 8: PARSE TABIX RESULTS
  # ==============================================================================
  
  cat("=== STEP 8: PARSING ANNOTATION DATA ===\n")
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
  # STEP 9: PROCESS ANNOTATION
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
  
  # Add chr prefix
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
  # STEP 10: COMPUTE OVERLAP
  # ==============================================================================
  
  cat("=== STEP 10: COMPUTING OVERLAP ===\n")
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
    cat("Created empty output file:", output_intersect, "\n")
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
  
  # Add annotation columns
  for (col in colnames(bedB)) {
    if (!col %in% c("Chromosome", "start", "end")) {
      intersect_df[[col]] <- bed2_hits[[col]]
    }
  }
  
  # Add highlight column
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
  cat("Total annotated variants:", nrow(intersect_df), "\n\n")
  
  # ==============================================================================
  # STEP 11: SAVE TSV OUTPUT
  # ==============================================================================
  
  cat("=== STEP 11: SAVING OUTPUT ===\n")
  
  write.table(intersect_df, output_intersect, 
              sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  cat("Saved intersect output:", output_intersect, "\n")
  
  # ==============================================================================
  # STEP 12: SAVE EXCEL OUTPUT WITH FORMATTING
  # ==============================================================================
  
  if (grepl("\\.xlsx?$", output_formatted, ignore.case = TRUE)) {
    cat("Creating formatted Excel file...\n")
    excel_time <- Sys.time()
    
    wb <- createWorkbook()
    addWorksheet(wb, "Annotated Variants")
    writeData(wb, "Annotated Variants", intersect_df)
    
    # Header style
    headerStyle <- createStyle(fontColour = "#FFFFFF", fgFill = "#4F81BD",
                               halign = "center", textDecoration = "Bold",
                               border = "TopBottomLeftRight")
    addStyle(wb, "Annotated Variants", headerStyle, rows = 1, 
             cols = 1:ncol(intersect_df), gridExpand = TRUE)

    # Set column widths
    setColWidths(wb, "Annotated Variants", cols = 1:ncol(intersect_df), widths = "auto")

    # Freeze header row
    freezePane(wb, "Annotated Variants", firstRow = TRUE)

    # Apply conditional formatting
      # ClinVar: red if not ".", green if "."
      if ("ClinVar" %in% colnames(intersect_df)) {
        clinvar_col <- which(colnames(intersect_df) == "ClinVar")
        green_style <- createStyle(fgFill = "#90EE90")
        red_style <- createStyle(fgFill = "#FFB6C1")
        for (i in 1:nrow(intersect_df)) {
          if (intersect_df$ClinVar[i] == ".") {
            addStyle(wb, "Annotated Variants", green_style, rows = i + 1, cols = clinvar_col)
          } else {
            addStyle(wb, "Annotated Variants", red_style, rows = i + 1, cols = clinvar_col)
          }
        }
      }
      # CADD_PHED: red if > 20, green otherwise
      if ("CADD_PHED" %in% colnames(intersect_df)) {
        cadd_col <- which(colnames(intersect_df) == "CADD_PHED")
        green_style <- createStyle(fgFill = "#90EE90")
        red_style <- createStyle(fgFill = "#FFB6C1")
        for (i in 1:nrow(intersect_df)) {
          val <- intersect_df$CADD_PHED[i]
          if (!is.na(val) && val > 20) {
            addStyle(wb, "Annotated Variants", red_style, rows = i + 1, cols = cadd_col)
          } else if (!is.na(val)) {
            addStyle(wb, "Annotated Variants", green_style, rows = i + 1, cols = cadd_col)
          }
        }
      }
      # MutationAssessor: red for H, yellow for M, green otherwise
      if ("MutationAssessor" %in% colnames(intersect_df)) {
        ma_col <- which(colnames(intersect_df) == "MutationAssessor")
        red_style <- createStyle(fgFill = "#FFB6C1")
        yellow_style <- createStyle(fgFill = "#FFFFE0")
        green_style <- createStyle(fgFill = "#90EE90")
        for (i in 1:nrow(intersect_df)) {
          val <- intersect_df$MutationAssessor[i]
          if (!is.na(val) && val == "H") {
            addStyle(wb, "Annotated Variants", red_style, rows = i + 1, cols = ma_col)
          } else if (!is.na(val) && val == "M") {
            addStyle(wb, "Annotated Variants", yellow_style, rows = i + 1, cols = ma_col)
          } else if (!is.na(val)) {
            addStyle(wb, "Annotated Variants", green_style, rows = i + 1, cols = ma_col)
          }
        }
      }
      # M_CAP: red for D, green otherwise
      if ("M_CAP" %in% colnames(intersect_df)) {
        mcap_col <- which(colnames(intersect_df) == "M_CAP")
        red_style <- createStyle(fgFill = "#FFB6C1")
        green_style <- createStyle(fgFill = "#90EE90")
        for (i in 1:nrow(intersect_df)) {
          val <- intersect_df$M_CAP[i]
          if (!is.na(val) && val == "D") {
            addStyle(wb, "Annotated Variants", red_style, rows = i + 1, cols = mcap_col)
          } else if (!is.na(val)) {
            addStyle(wb, "Annotated Variants", green_style, rows = i + 1, cols = mcap_col)
          }
        }
      }
      # AF_gnomAD: red if < 0.5, green otherwise
      if ("AF_gnomAD" %in% colnames(intersect_df)) {
        af_col <- which(colnames(intersect_df) == "AF_gnomAD")
        red_style <- createStyle(fgFill = "#FFB6C1")
        green_style <- createStyle(fgFill = "#90EE90")
        for (i in 1:nrow(intersect_df)) {
          val <- intersect_df$AF_gnomAD[i]
          if (!is.na(val) && val < 0.5) {
            addStyle(wb, "Annotated Variants", red_style, rows = i + 1, cols = af_col)
          } else if (!is.na(val)) {
            addStyle(wb, "Annotated Variants", green_style, rows = i + 1, cols = af_col)
          }
        }
      }
      # Highlight important variants (start/end columns)
      if ("highlight_important" %in% colnames(intersect_df)) {
        start_col <- which(colnames(intersect_df) == "start")
        end_col <- which(colnames(intersect_df) == "end")
        yellow_style <- createStyle(fgFill = "#FFFF00", fontColour = "#FF0000")
        white_style <- createStyle(fgFill = "#FFFFFF", fontColour = "#00FF00")
        for (i in 1:nrow(intersect_df)) {
          if (intersect_df$highlight_important[i]) {
            addStyle(wb, "Annotated Variants", yellow_style, rows = i + 1, cols = start_col)
            addStyle(wb, "Annotated Variants", yellow_style, rows = i + 1, cols = end_col)
          } else {
            addStyle(wb, "Annotated Variants", white_style, rows = i + 1, cols = start_col)
            addStyle(wb, "Annotated Variants", white_style, rows = i + 1, cols = end_col)
          }
        }
      }
    
    saveWorkbook(wb, output_formatted, overwrite = TRUE)
    cat("Saved formatted output:", output_formatted, "\n")

    cat(paste("Excel formatting took:", 
              round(difftime(Sys.time(), excel_time, units = "secs"), 2), 
              "seconds\n"))
  }

  cat("\n=== DONE ===\n")

  invisible(intersect_df)
}