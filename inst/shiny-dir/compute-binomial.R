#============================================================================
# compute-binomial.R - BINOMIAL DISTRIBUTION ANALYSIS (INTEGRATED WITH VARIANTS)
# ============================================================================

# ============================================================================
# REACTIVE: Get coverage at specific genomic position FROM VARIANTS
# ============================================================================

df_subset <- reactive({
  cat("\n=== BINOMIAL: df_subset called ===\n")
  
  # Validazione input
  if (is.null(input$start_gp) || input$start_gp == "") {
    cat("ERROR: No genomic position provided\n")
    return(NULL)
  }

  cat("Input genomic position:", input$start_gp, "\n")
  
  # CAMBIAMENTO CHIAVE: Usa i dati delle varianti invece di mysample()
  variant_data <- tryCatch({
    uncover_maxaf_data()  # Oppure annotated_variants_data() se preferisci
  }, error = function(e) {
    cat("ERROR getting variant data:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(variant_data) || nrow(variant_data) == 0) {
    cat("ERROR: No variant data available\n")
    return(NULL)
  }
  
  cat("Variant data rows:", nrow(variant_data), "\n")
  cat("Variant data columns:", paste(colnames(variant_data), collapse=", "), "\n")
  
  # Converti input a numeric
  genomic_pos <- suppressWarnings(as.numeric(input$start_gp))
  
  if (is.na(genomic_pos)) {
    cat("ERROR: Cannot convert", input$start_gp, "to numeric\n")
    return(NULL)
  }
  
  cat("Looking for position:", genomic_pos, "\n")
  
  # Cerca la posizione esatta o la più vicina
  if (!genomic_pos %in% variant_data$start) {
    cat("WARNING: Exact position", genomic_pos, "not found\n")
    
    # Trova la posizione più vicina
    closest_idx <- which.min(abs(variant_data$start - genomic_pos))
    closest_pos <- variant_data$start[closest_idx]
    distance <- abs(closest_pos - genomic_pos)
    
    cat("Using closest variant at position:", closest_pos, "\n")
    cat("Distance:", distance, "bp\n")
    
    # Se troppo distante, avvisa
    if (distance > 1000) {
      cat("WARNING: Closest variant is", distance, "bp away!\n")
    }
    
    genomic_pos <- closest_pos
  }
  
  # Subset data
  position_data <- variant_data %>%
    dplyr::filter(start == genomic_pos)
  
  cat("Filtered rows:", nrow(position_data), "\n")
  
  if (nrow(position_data) == 0) {
    cat("ERROR: No data after filtering\n")
    return(NULL)
  }
  
  # Extract coverage e altre info
  cov_value <- as.numeric(position_data$coverage[1])
  af_gnomad <- as.numeric(position_data$AF_gnomAD[1])
  gene_name <- as.character(position_data$GENENAME[1])
  ref_allele <- as.character(position_data$REF[1])
  alt_allele <- as.character(position_data$ALT[1])
  chromosome <- as.character(position_data$seqnames[1])
  
  cat("Coverage value:", cov_value, "\n")
  cat("AF_gnomAD:", af_gnomad, "\n")
  cat("Gene:", gene_name, "\n")
  cat("Variant:", ref_allele, ">", alt_allele, "\n")
  
  if (is.na(cov_value) || cov_value <= 0) {
    cat("ERROR: Invalid coverage\n")
    return(NULL)
  }
  
  cat("SUCCESS: Returning data for position", genomic_pos, "\n")
  
  # Ritorna lista completa di informazioni
  return(list(
    coverage = cov_value,
    position_used = genomic_pos,
    position_requested = as.numeric(input$start_gp),
    af_gnomad = af_gnomad,
    gene = gene_name,
    ref = ref_allele,
    alt = alt_allele,
    chromosome = chromosome
  ))
})

# ============================================================================
# OUTPUT: Binomial Distribution Barplot
# ============================================================================

output$bd <- renderPlot({
  cat("\n=== RENDERING BD PLOT ===\n")
  
  result <- df_subset()
  
  if (is.null(result)) {
    plot.new()
    text(0.5, 0.5, "No variant data available\nPlease enter a valid genomic position", 
         cex = 1.5, col = "red")
    return()
  }
  
  coverage <- result$coverage
  pos_used <- result$position_used
  pos_requested <- result$position_requested
  gene_name <- result$gene
  
  cat("Coverage received:", coverage, "\n")
  cat("Coverage class:", class(coverage), "\n")
  
  if (is.na(coverage) || coverage <= 0) {
    plot.new()
    text(0.5, 0.5, paste0("Invalid coverage: ", coverage), 
         cex = 1.5, col = "red")
    return()
  }
  
  if (is.null(input$p) || is.na(input$p) || input$p <= 0 || input$p > 1) {
    plot.new()
    text(0.5, 0.5, "Invalid allele fraction\nMust be between 0 and 1", 
         cex = 1.5, col = "red")
    return()
  }
  
  cat("Calculating CI with size =", coverage, "prob =", input$p, "\n")
  
  # Calculate CI
  ci <- tryCatch({
    qbinom(p = c(0.025, 0.975), size = as.integer(coverage), prob = as.numeric(input$p))
  }, error = function(e) {
    cat("ERROR in qbinom:", e$message, "\n")
    return(c(0, coverage))
  })
  
  cat("CI:", ci, "\n")
  
  # Calculate probabilities
  x_values <- 0:coverage
  probabilities <- dbinom(x = x_values, size = as.integer(coverage), prob = as.numeric(input$p))
  
  # Nota se posizione diversa
  pos_note <- if (pos_used != pos_requested) {
    paste0("\n(Using closest variant at ", pos_used, ")")
  } else {
    ""
  }
  
  # Informazioni variante
  variant_info <- paste0(
    "Gene: ", ifelse(is.na(gene_name) || gene_name == ".", "Unknown", gene_name),
    " | ", result$ref, ">", result$alt
  )
  
  # Create barplot
  barplot(
    probabilities,
    names.arg = x_values,
    main = paste0("Binomial Distribution\n",
                  "Position: ", pos_requested, pos_note, "\n",
                  variant_info, "\n",
                  "Coverage: ", coverage, 
                  " | Expected AF: ", input$p,
                  "\n95% CI: [", ci[1], ", ", ci[2], "]"),
    xlab = "Number of Variant Reads",
    ylab = "Probability",
    col = ifelse(x_values >= input$num_all, "lightcoral", "lightblue"),
    border = "white",
    las = 1,
    cex.main = 0.85
  )
  
  abline(v = input$num_all, col = "red", lwd = 2, lty = 2)
  abline(v = ci[1], col = "blue", lwd = 1, lty = 3)
  abline(v = ci[2], col = "blue", lwd = 1, lty = 3)
  
  legend("topright", 
         legend = c(paste0("Required reads (", input$num_all, ")"), "95% CI bounds"),
         col = c("red", "blue"),
         lty = c(2, 3),
         lwd = c(2, 1),
         bty = "n")
})

# ============================================================================
# OUTPUT: Cumulative Distribution Function Plot
# ============================================================================

output$pbinom <- renderPlot({
  cat("\n=== RENDERING PBINOM PLOT ===\n")
  
  result <- df_subset()
  
  if (is.null(result)) {
    plot.new()
    text(0.5, 0.5, "No valid variant data", cex = 1.5, col = "red")
    return()
  }
  
  coverage <- result$coverage
  
  if (is.na(coverage) || coverage <= 0) {
    plot.new()
    text(0.5, 0.5, "Invalid coverage", cex = 1.5, col = "red")
    return()
  }
  
  if (is.null(input$p) || is.na(input$p) || input$p <= 0 || input$p > 1) {
    plot.new()
    text(0.5, 0.5, "Invalid allele fraction", cex = 1.5, col = "red")
    return()
  }
  
  if (is.null(input$num_all) || is.na(input$num_all) || input$num_all < 0) {
    plot.new()
    text(0.5, 0.5, "Invalid variant reads threshold", cex = 1.5, col = "red")
    return()
  }
  
  x_values <- 0:coverage
  cum_probs <- pbinom(q = x_values, size = as.integer(coverage), prob = as.numeric(input$p))
  
  prob_more <- 1 - pbinom(q = as.integer(input$num_all), 
                          size = as.integer(coverage), 
                          prob = as.numeric(input$p))
  
  # Info variante per il titolo
  gene_name <- ifelse(is.na(result$gene) || result$gene == ".", "Unknown", result$gene)
  variant_info <- paste0(gene_name, " (", result$ref, ">", result$alt, ")")
  
  plot(x_values, cum_probs, 
       type = 'h', 
       lwd = 2,
       col = "steelblue",
       xlab = "Number of Variant Reads",
       ylab = "Cumulative Probability P(X ≤ x)",
       main = paste0("Cumulative Distribution Function\n",
                     variant_info, "\n",
                     "P(X > ", input$num_all, ") = ", 
                     round(prob_more * 100, 2), "%"),
       las = 1)
  
  abline(v = input$num_all, col = "red", lwd = 2, lty = 2)
  
  if (input$num_all <= coverage) {
    prob_at <- pbinom(q = as.integer(input$num_all), 
                      size = as.integer(coverage), 
                      prob = as.numeric(input$p))
    abline(h = prob_at, col = "red", lwd = 1, lty = 3)
    
    text(x = coverage * 0.7, 
         y = prob_at + 0.05,
         labels = paste0("P(X ≤ ", input$num_all, ") = ", round(prob_at * 100, 2), "%"),
         col = "red",
         pos = 3)
  }
  
  grid(col = "gray90")
})

# ============================================================================
# OUTPUT: Interpretation Text
# ============================================================================

output$ci <- renderText({
  cat("\n=== RENDERING CI TEXT ===\n")
  
  result <- df_subset()
  
  if (is.null(result)) {
    return("<span style='color:red'><strong>Please enter a valid genomic position</strong></span>")
  }
  
  coverage <- result$coverage
  
  if (is.na(coverage) || coverage <= 0) {
    return("<span style='color:red'><strong>Invalid coverage data</strong></span>")
  }
  
  if (is.null(input$p) || is.na(input$p) || input$p <= 0 || input$p > 1) {
    return("<span style='color:red'><strong>Invalid allele fraction (must be 0-1)</strong></span>")
  }
  
  if (is.null(input$num_all) || is.na(input$num_all)) {
    return("<span style='color:red'><strong>Please enter variant reads threshold</strong></span>")
  }
  
  if (input$num_all > coverage) {
    return(paste0("<span style='color:red'><strong>ERROR:</strong> Required reads (", 
                  input$num_all, ") exceeds coverage (", coverage, ")</span>"))
  }
  
  ci <- qbinom(p = c(0.025, 0.975), size = as.integer(coverage), prob = as.numeric(input$p))
  prob_more <- 1 - pbinom(q = as.integer(input$num_all) - 1, 
                          size = as.integer(coverage), 
                          prob = as.numeric(input$p))
  
  # Info variante
  gene_display <- ifelse(is.na(result$gene) || result$gene == ".", "Unknown gene", result$gene)
  variant_info <- paste0(
    "<strong>Variant:</strong> ", 
    gene_display, " ", 
    result$ref, ">", result$alt, 
    " @ chr", result$chromosome, ":", result$position_used, "<br>"
  )
  
  # Info gnomAD
  gnomad_info <- if (!is.na(result$af_gnomad)) {
    paste0("<strong>gnomAD AF:</strong> ", signif(result$af_gnomad, 3), "<br>")
  } else {
    "<strong>gnomAD AF:</strong> Not available<br>"
  }
  
  # Nota se posizione diversa da quella richiesta
  position_note <- if (result$position_used != result$position_requested) {
    paste0("<strong>Note:</strong> Using closest variant (requested: ", 
           result$position_requested, 
           ", distance: ", 
           abs(result$position_used - result$position_requested), 
           " bp)<br>")
  } else {
    ""
  }
  
  if (input$num_all < ci[1]) {
    color <- "red"
    interpretation <- paste0(
      "⚠️ The required ", input$num_all, " variant reads is BELOW the 95% CI [",
      ci[1], ", ", ci[2], "]. ",
      "<strong>This suggests you might be MISSING true variants!</strong>"
    )
  } else if (input$num_all > ci[2]) {
    color <- "red"
    interpretation <- paste0(
      "⚠️ The required ", input$num_all, " variant reads is ABOVE the 95% CI [",
      ci[1], ", ", ci[2], "]. ",
      "<strong>This threshold is too stringent.</strong>"
    )
  } else {
    color <- "green"
    interpretation <- paste0(
      "✅ The required ", input$num_all, " variant reads falls WITHIN the 95% CI [",
      ci[1], ", ", ci[2], "]. ",
      "<strong>This is a reasonable threshold.</strong>"
    )
  }
  
  paste0("<div style='color:", color, "; padding: 10px; background-color: #f9f9f9; border-left: 4px solid ", color, ";'>",
         variant_info,
         gnomad_info,
         position_note,
         "<strong>Coverage:</strong> ", coverage, " reads<br>",
         "<strong>Expected AF:</strong> ", input$p * 100, "%<br>",
         "<strong>95% CI:</strong> [", ci[1], ", ", ci[2], "] reads<br>",
         "<strong>P(X ≥ ", input$num_all, "):</strong> ", round(prob_more * 100, 2), "%<br><br>",
         "<strong>Interpretation:</strong><br>",
         interpretation,
         "</div>")
})


