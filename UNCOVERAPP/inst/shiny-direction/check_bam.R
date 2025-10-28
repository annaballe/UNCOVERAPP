#!/usr/bin/env Rscript

# BAM File Integrity Checker
# This script comprehensively checks BAM file integrity and provides detailed diagnostics

library(Rsamtools)
library(GenomicAlignments)
library(dplyr)
library(tibble)

# Function to check BAM file integrity
check_bam_integrity <- function(file_path, verbose = TRUE) {
  
  if (verbose) cat("Checking:", basename(file_path), "\n")
  
  result <- list(
    file = basename(file_path),
    full_path = file_path,
    exists = FALSE,
    is_bam = FALSE,
    can_open = FALSE,
    has_header = FALSE,
    has_sq_lines = FALSE,
    has_reads = FALSE,
    is_sorted = FALSE,
    is_indexed = FALSE,
    file_size_mb = 0,
    error_message = "",
    read_count = 0,
    chromosome_count = 0,
    warnings = character(0)
  )
  
  # Check if file exists
  if (!file.exists(file_path)) {
    result$error_message <- "File does not exist"
    return(result)
  }
  
  result$exists <- TRUE
  result$file_size_mb <- round(file.size(file_path) / 1024^2, 2)
  
  # Check if it's actually a BAM file by extension and magic bytes
  if (!grepl("\\.bam$", file_path, ignore.case = TRUE)) {
    result$error_message <- "Not a BAM file (wrong extension)"
    return(result)
  }
  
  # Check BAM magic bytes
  tryCatch({
    con <- file(file_path, "rb")
    magic_bytes <- readBin(con, "raw", 4)
    close(con)
    
    # BAM files start with specific magic bytes
    if (!identical(magic_bytes, as.raw(c(0x42, 0x41, 0x4d, 0x01)))) {
      result$error_message <- "Not a valid BAM file (wrong magic bytes)"
      return(result)
    }
  }, error = function(e) {
    result$error_message <- paste("Cannot read file:", e$message)
    return(result)
  })
  
  result$is_bam <- TRUE
  
  # Try to open with Rsamtools
  tryCatch({
    bf <- BamFile(file_path)
    open(bf)
    result$can_open <- TRUE
    
    # Check header
    header <- scanBamHeader(bf)
    if (length(header) > 0) {
      result$has_header <- TRUE
      
      # Check for @SQ lines (sequence dictionary)
      if ("targets" %in% names(header[[1]]) && length(header[[1]]$targets) > 0) {
        result$has_sq_lines <- TRUE
        result$chromosome_count <- length(header[[1]]$targets)
      } else {
        result$warnings <- c(result$warnings, "No @SQ lines in header")
      }
    }
    
    # Check if file has reads
    reads <- scanBam(bf, param = ScanBamParam(what = "qname"))
    if (length(reads[[1]]$qname) > 0) {
      result$has_reads <- TRUE
      result$read_count <- length(reads[[1]]$qname)
    }
    
    # Check if sorted
    if (result$has_header) {
      so_tag <- header[[1]]$text[grepl("^@HD", header[[1]]$text)]
      if (length(so_tag) > 0) {
        if (grepl("SO:coordinate", so_tag)) {
          result$is_sorted <- TRUE
        } else if (grepl("SO:queryname", so_tag)) {
          result$warnings <- c(result$warnings, "Sorted by queryname, not coordinate")
        }
      }
    }
    
    close(bf)
    
  }, error = function(e) {
    result$error_message <- paste("Cannot open BAM file:", e$message)
    return(result)
  })
  
  # Check if indexed
  index_file <- paste0(file_path, ".bai")
  if (file.exists(index_file)) {
    result$is_indexed <- TRUE
  } else {
    result$warnings <- c(result$warnings, "No index file (.bai) found")
  }
  
  return(result)
}

# Function to check multiple BAM files
check_multiple_bam_files <- function(file_paths, verbose = TRUE) {
  
  cat("=== BAM FILE INTEGRITY CHECK ===\n")
  cat("Checking", length(file_paths), "files...\n\n")
  
  results <- list()
  
  for (i in seq_along(file_paths)) {
    if (verbose) cat(paste0("[", i, "/", length(file_paths), "] "))
    results[[i]] <- check_bam_integrity(file_paths[i], verbose = verbose)
    if (verbose) cat("\n")
  }
  
  # Convert to data frame for easier viewing
  df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      file = x$file,
      exists = x$exists,
      is_bam = x$is_bam,
      can_open = x$can_open,
      has_header = x$has_header,
      has_sq_lines = x$has_sq_lines,
      has_reads = x$has_reads,
      is_sorted = x$is_sorted,
      is_indexed = x$is_indexed,
      size_mb = x$file_size_mb,
      read_count = x$read_count,
      chrom_count = x$chromosome_count,
      error = x$error_message,
      warnings = paste(x$warnings, collapse = "; "),
      stringsAsFactors = FALSE
    )
  }))
  
  return(list(results = results, summary = df))
}

# Function to get file list from directory
get_potential_bam_files <- function(directory = ".", pattern = NULL) {
  
  if (is.null(pattern)) {
    # Get all files that might be BAM files
    all_files <- list.files(directory, full.names = TRUE, recursive = FALSE)
    
    # Filter for potential BAM files (including those without proper extensions)
    potential_bams <- all_files[file.info(all_files)$size > 100]  # At least 100 bytes
    
    return(potential_bams)
  } else {
    return(list.files(directory, pattern = pattern, full.names = TRUE))
  }
}

# Main execution function
run_bam_integrity_check <- function(directory = ".", pattern = "\\.bam$") {
  
  cat("Scanning directory:", directory, "\n")
  
  # Get file list
  if (is.null(pattern)) {
    files <- get_potential_bam_files(directory)
  } else {
    files <- list.files(directory, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  }
  
  if (length(files) == 0) {
    cat("No files found matching pattern!\n")
    return(NULL)
  }
  
  cat("Found", length(files), "potential BAM files:\n")
  cat(paste("-", basename(files)), sep = "\n")
  cat("\n")
  
  # Check integrity
  check_results <- check_multiple_bam_files(files)
  
  # Print summary
  cat("\n=== SUMMARY ===\n")
  print(check_results$summary)
  
  # Print detailed results for failed files
  failed_files <- check_results$results[sapply(check_results$results, function(x) !x$can_open)]
  if (length(failed_files) > 0) {
    cat("\n=== FAILED FILES DETAILS ===\n")
    for (failed in failed_files) {
      cat("\nFile:", failed$file, "\n")
      cat("Error:", failed$error_message, "\n")
      cat("Size:", failed$file_size_mb, "MB\n")
      if (length(failed$warnings) > 0) {
        cat("Warnings:", paste(failed$warnings, collapse = ", "), "\n")
      }
    }
  }
  
  # Print recommendations
  valid_files <- check_results$results[sapply(check_results$results, function(x) x$can_open)]
  cat("\n=== RECOMMENDATIONS ===\n")
  cat("Valid BAM files found:", length(valid_files), "\n")
  cat("Invalid BAM files found:", length(failed_files), "\n")
  
  if (length(valid_files) > 0) {
    cat("\nValid BAM files:\n")
    for (valid in valid_files) {
      cat("-", valid$file, 
          paste0("(", valid$read_count, " reads, ", valid$chromosome_count, " chromosomes)"), "\n")
    }
  }
  
  return(check_results)
}

# Example usage:
# Run this to check all potential BAM files in current directory
cat("To run the integrity check, use:\n")
cat("results <- run_bam_integrity_check()\n")
cat("Or for a specific directory:\n")
cat("results <- run_bam_integrity_check('/path/to/your/bam/files')\n")
cat("\nTo check only .bam files:\n")
cat("results <- run_bam_integrity_check(pattern = '\\.bam$')\n")
cat("\nTo check ALL files (including those without .bam extension):\n")
cat("results <- run_bam_integrity_check(pattern = NULL)\n")