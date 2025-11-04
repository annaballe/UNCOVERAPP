# ============================================================================
# compute-reactiveDF.R - GESTIONE DATI BASE
# ============================================================================

# Database reference genome
txdb <- reactive({
  if (input$UCSC_Genome == "hg19") {
    TxDb.Hsapiens.UCSC.hg19.knownGene
  } else {
    TxDb.Hsapiens.UCSC.hg38.knownGene
  }
})

# ============================================================================
# DATA SOURCE MANAGEMENT
# ============================================================================

data_source <- reactiveVal("none")
raw_upload <- reactiveVal(NULL)

observeEvent(input$file1, {
  req(input$file1)
  raw_upload(input$file1)
  data_source("manual")
})

observeEvent(input$pileup, {
  req(coverage_input())
  data_source("pileup")
})

# ============================================================================
# MAIN DATA REACTIVE
# ============================================================================

mydata <- reactive({
  source <- data_source()
  
  if (source == "none") return(NULL)
  
  if (source == "manual") {
    req(raw_upload())
    file_info <- raw_upload()
    
    cat("\n=== PROCESSING MANUAL UPLOAD ===\n")
    cat("File:", file_info$name, "\n")
    
    tmp <- read.table(file_info$datapath,
                      header = input$header, 
                      stringsAsFactors = FALSE)
    
    cat("Raw dimensions:", dim(tmp), "\n")
    
    if (ncol(tmp) >= 3) {
      colnames(tmp)[1:3] <- c("chromosome", "start", "end")
    }
    
    tmp$chromosome <- paste0("chr", sub("^chr", "", as.character(tmp$chromosome)))
    tmp$start <- as.integer(tmp$start)
    tmp$end <- as.integer(tmp$end)
    
    if (any(grepl("^sample_", names(tmp)))) {
      cat("Found existing sample_* columns\n")
      return(tmp)
    }
    
    if (ncol(tmp) >= 4) {
      sample_label <- tools::file_path_sans_ext(basename(file_info$name))
      cov <- suppressWarnings(as.integer(tmp[[4]]))
      
      new_df <- data.frame(
        chromosome = tmp$chromosome,
        start = tmp$start,
        end = tmp$end,
        cov = cov,
        stringsAsFactors = FALSE
      )
      colnames(new_df)[4] <- paste0("sample_", sample_label)
      
      cat("Created sample column:", colnames(new_df)[4], "\n")
      return(new_df)
    }
    
    return(NULL)
  }
  
  if (source == "pileup") {
    req(coverage_input())
    
    cat("\n=== PROCESSING PILEUP DATA ===\n")
    
    tmp_pileup <- coverage_input()
    cat("Coverage dimensions:", dim(tmp_pileup), "\n")
    
    colnames(tmp_pileup)[1:3] <- c("chromosome", "start", "end")
    tmp_pileup$chromosome <- paste0("chr", sub("^chr", "", as.character(tmp_pileup$chromosome)))
    tmp_pileup$chromosome <- as.character(tmp_pileup$chromosome)
    tmp_pileup$start <- as.integer(tmp_pileup$start)
    tmp_pileup$end <- as.integer(tmp_pileup$end)
    
    data_cols <- setdiff(colnames(tmp_pileup), c("chromosome", "start", "end"))
    ncols <- length(data_cols)
    samples <- name_sample()
    k <- length(samples)
    
    cat("Data columns:", ncols, "| Samples:", k, "\n")
    
    if (ncols == k && k > 0) {
      colnames(tmp_pileup)[-1:-3] <- paste0("sample_", samples)
      cat("Applied BED naming\n")
    } else if (ncols == 2 * k && k > 0) {
      interleaved <- as.vector(rbind(
        paste0("sample_", samples), 
        paste0("nucleotide_", samples)
      ))
      colnames(tmp_pileup)[-1:-3] <- interleaved
      cat("Applied BAM naming\n")
    } else {
      colnames(tmp_pileup)[-1:-3] <- paste0("sample_", seq_len(ncols))
      cat("WARNING: Using fallback naming\n")
    }
    
    cat("Final columns:", paste(colnames(tmp_pileup), collapse=", "), "\n")
    
    cov_cols <- grep("^sample_", colnames(tmp_pileup), value = TRUE)
    tmp_pileup[cov_cols] <- lapply(tmp_pileup[cov_cols], as.integer)
    
    return(tmp_pileup)
  }
  
  return(NULL)
})

# ============================================================================
# SAMPLE EXTRACTION
# ============================================================================

mysample <- reactive({
  dat <- mydata()
  req(dat, input$Sample)
  
  sel <- input$Sample
  if (sel == "") return(NULL)
  
  cat("\n=== EXTRACTING SAMPLE ===\n")
  cat("Requested:", sel, "\n")
  
  sel_base <- tools::file_path_sans_ext(sel)
  sel_clean <- sub("^sample_", "", sel)
  sel_clean_base <- tools::file_path_sans_ext(sel_clean)
  
  candidates <- unique(c(
    sel, sel_base, sel_clean, sel_clean_base,
    paste0("sample_", sel),
    paste0("sample_", sel_base),
    paste0("sample_", sel_clean),
    paste0("sample_", sel_clean_base)
  ))
  
  col <- intersect(candidates, names(dat))
  
  if (length(col) == 0) {
    cat("ERROR: No match found\n")
    cat("Available:", paste(names(dat), collapse=", "), "\n")
    return(NULL)
  }
  
  col <- col[1]
  cat("Selected:", col, "\n")
  
  out <- dat[, c("chromosome", "start", "end", col), drop = FALSE]
  names(out)[4] <- "coverage"
  
  cat("Extracted", nrow(out), "rows\n")
  
  return(out)
})

# ============================================================================
# FILTERING REACTIVES
# ============================================================================

filtered_low <- reactive({
  req(mysample(), input$coverage_co, input$Chromosome)
  
  df <- mysample()
  thr <- input$coverage_co
  chr_val <- input$Chromosome
  
  cat("\n=== FILTERING LOW COVERAGE ===\n")
  cat("Chr:", chr_val, "| Threshold:", thr, "\n")
  
  result <- if (identical(thr, "all")) {
    dplyr::filter(df, chromosome == chr_val)
  } else {
    dplyr::filter(df, chromosome == chr_val, coverage <= as.numeric(thr))
  }
  
  cat("Result:", nrow(result), "rows\n")
  return(result)
})

filtered_high <- reactive({
  req(mysample(), input$coverage_co, input$Chromosome)
  
  df <- mysample()
  thr <- input$coverage_co
  chr_val <- input$Chromosome
  
  result <- if (identical(thr, "all")) {
    dplyr::filter(df, chromosome == chr_val)
  } else {
    dplyr::filter(df, chromosome == chr_val, coverage > as.numeric(thr))
  }
  
  return(result)
})

# ALIAS per compatibilità con compute-annotation.R
filtered_low_nucl <- filtered_low