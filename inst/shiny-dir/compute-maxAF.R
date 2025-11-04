# ============================================================================
# compute-maxAF.R - MAX ALLELE FREQUENCY CALCULATOR
# ============================================================================

# ============================================================================
# CALCOLO maxAF e maxAC
# ============================================================================

data <- reactive({
  myPrev <- 1 / input$prev
  
  if (input$inh == "monoallelic") {
    myMaxAF <- (1/2) * myPrev * input$hetA * input$hetG * (1/input$pen)
  } else if (input$inh == "biallelic") {
    myMaxAF <- sqrt(myPrev) * input$hetA * sqrt(input$hetG) * (1/sqrt(input$pen))
  } else {
    myMaxAF <- NA
  }
  
  myMaxAC <- qpois(
    p = as.numeric(input$CI),
    lambda = (input$popSize) * (myMaxAF)
  )
  
  return(list(myMaxAF, myMaxAC))
})

# ============================================================================
# DATAFRAME PULITO: Filtra per maxAF threshold
# ============================================================================

uncover_maxaf_data <- reactive({
  cat("\n=== MAX AF FILTERING ===\n")
  
  # Prendi i dati annotati base
  base_data <- annotated_variants_data()
  
  if (is.null(base_data) || nrow(base_data) == 0) {
    return(NULL)
  }
  
  # Calcola maxAF threshold
  maxAF_threshold <- signif(data()[[1]], 3)
  cat("MaxAF threshold:", maxAF_threshold, "\n")
  
  # Filtra per maxAF
  filtered <- base_data %>%
    dplyr::filter(
      !is.na(AF_gnomAD),
      AF_gnomAD < maxAF_threshold
    )
  
  cat("Variants under maxAF:", nrow(filtered), "\n")
  
  # CRITICAL FIX: Se 0 righe, ritorna dataframe vuoto subito
  if (nrow(filtered) == 0) {
    cat("WARNING: No variants pass maxAF filter (all AF >= ", maxAF_threshold, ")\n", sep="")
    # Ritorna dataframe vuoto con struttura corretta
    return(data.frame(
      seqnames = character(),
      start = integer(),
      end = integer(),
      coverage = numeric(),
      REF = character(),
      ALT = character(),
      dbsnp = character(),
      GENENAME = character(),
      PROTEIN_ensembl = character(),
      MutationAssessor = character(),
      SIFT = character(),
      Polyphen2 = character(),
      M_CAP = character(),
      CADD_PHED = numeric(),
      AF_gnomAD = numeric(),
      ClinVar = character(),
      clinvar_MedGen_id = character(),
      HGVSc_VEP = character(),
      HGVSp_VEP = character(),
      highlight_important = logical(),
      highlight_maxaf = logical(),
      counts = numeric(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Se arriviamo qui, abbiamo dati: aggiungi colonne helper
  filtered$highlight_maxaf <- grepl("H|M", filtered$MutationAssessor) & 
                              filtered$ClinVar != "." & 
                              filtered$M_CAP != "TRUE" &
                              !is.na(filtered$AF_gnomAD) & 
                              filtered$AF_gnomAD < maxAF_threshold
  
  # Aggiungi counts column se non esiste
  if (!"counts" %in% colnames(filtered)) {
    filtered$counts <- NA
  }
  
  cat("Returning", nrow(filtered), "variants\n")
  return(filtered)
})

# ============================================================================
# TABELLA FORMATTATA CON DT (invece di condformat)
# ============================================================================

uncover_maxaf <- reactive({
  cat("\n=== MAX AF TABLE (DT WIDGET) ===\n")
  
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Filtering by maxAF...", value = 0.5)
  
  data <- uncover_maxaf_data()
  
  if (is.null(data) || nrow(data) == 0) {
    cat("WARNING: No variants pass maxAF filter!\n")
    # Restituisci un messaggio informativo invece di NULL
    return(
      DT::datatable(
        data.frame(Message = "No variants found with AF < maxAF threshold"),
        options = list(dom = 't'),
        rownames = FALSE
      )
    )
  }
  
  cat("Formatting", nrow(data), "variants with maxAF filter\n")
  
  # Calcola maxAF per mostrarlo nella tabella
  maxAF_value <- signif(data()[[1]], 3)
  
  dt <- DT::datatable(
    data,
    options = list(
      pageLength = 25,
      scrollX = TRUE,
      scrollY = "600px",
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel'),
      columnDefs = list(
        list(targets = which(colnames(data) == "highlight_maxaf") - 1, visible = FALSE),
        list(targets = which(colnames(data) == "highlight_important") - 1, visible = FALSE)
      )
    ),
    extensions = 'Buttons',
    rownames = FALSE,
    caption = htmltools::tags$caption(
      style = 'caption-side: top; text-align: center; color: black; font-size: 150%;',
      paste0('Variants with AF < maxAF (', maxAF_value, ')')
    )
  ) %>%
    # ClinVar: rosso se presente, verde se assente
    DT::formatStyle(
      'ClinVar',
      backgroundColor = DT::styleEqual(
        c('.'), 
        c('lightgreen'),
        default = 'lightcoral'
      )
    ) %>%
    # MutationAssessor: rosso per H, giallo per M, verde altri
    DT::formatStyle(
      'MutationAssessor',
      backgroundColor = DT::styleEqual(
        c('H', 'M'), 
        c('lightcoral', 'lightyellow'),
        default = 'lightgreen'
      )
    ) %>%
    # M_CAP: rosso se deleterious
    DT::formatStyle(
      'M_CAP',
      backgroundColor = DT::styleEqual(
        c('D'), 
        c('lightcoral'),
        default = 'lightgreen'
      )
    ) %>%
    # AF_gnomAD: gradient based on maxAF
    DT::formatStyle(
      'AF_gnomAD',
      background = DT::styleColorBar(range(c(0, maxAF_value), na.rm = TRUE), 'lightblue'),
      backgroundSize = '100% 90%',
      backgroundRepeat = 'no-repeat',
      backgroundPosition = 'center'
    ) %>%
    # Highlight important variants (yellow background)
    DT::formatStyle(
      c('start', 'end'),
      'highlight_maxaf',
      backgroundColor = DT::styleEqual(
        c(TRUE, FALSE),
        c('yellow', 'white')
      )
    ) %>%
    # Color text (red if important, green otherwise)
    DT::formatStyle(
      c('start', 'end'),
      'highlight_maxaf',
      color = DT::styleEqual(
        c(TRUE, FALSE),
        c('red', 'green')
      )
    )
  
  cat("maxAF DT widget created\n")
  return(dt)
})