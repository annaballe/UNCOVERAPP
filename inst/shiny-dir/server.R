server <- function (input, output, session){
  suppressPackageStartupMessages({
    require(dplyr)
    require(Gviz)
    require(Homo.sapiens)
    require(stringr)
    require(shinyjs)
    require(bedr)
    require(EnsDb.Hsapiens.v75)
    require(EnsDb.Hsapiens.v86)
    require(DT)
    require(shinycssloaders)
  })

  options(shiny.maxRequestSize=30*1024^2)

  script1 <- system.file(
    "extdata",
    "Rpreprocessing.R",
    package = "uncoverappLib")

  # ============================================================================
  # DOWNLOAD DEPENDENCY SCRIPT
  # ============================================================================
  
  output$dependence = downloadHandler(
    filename="Rpreprocessing.R",
    content=function(file){
      file.copy(script1,file)
    }
  )
  
  # ============================================================================
  # SOURCE REACTIVE SCRIPTS
  # ============================================================================
  
  source('compute-preprocess.R', local= TRUE)
  source('compute-reactiveDF.R', local= TRUE)      # mydata(), mysample(), filtered_low()
  source('compute-annotation.R', local=TRUE)       # intBED(), annotated_variants_data(), condform_table()
  source('compute-tables.R', local= TRUE)
  source('compute-plots.R', local=TRUE)
  source('compute-maxAF.R', local=TRUE)            # data(), uncover_maxaf_data(), uncover_maxaf()
  source('compute-binomial.R', local=TRUE)
  
  # ============================================================================
  # OUTPUT: COVERAGE INPUT TABLE
  # ============================================================================
  
  output$input1 <- renderDataTable({
    options(shiny.sanitize.errors = TRUE)
      start_time <- Sys.time()   # Start time
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "table construction in progress",
                 detail = 'This may take a while', value = 0)
    Sys.sleep(0.1)
    validate(
      need(try(!is.null(coverage_input())), "please upload a file with HGNC
      gene names and absolute path(s) to BAM file"))
  
    df <- coverage_input()
  
    end_time <- Sys.time()     # ⏱️ end
    print(paste("Tempo di esecuzione:", round(end_time - start_time, 2)))  # console
  
    df
  })

  # ============================================================================
  # DOWNLOAD: STATISTICAL SUMMARY
  # ============================================================================
  
  output$summary <- downloadHandler(
    filename = function() {
      paste('statistical_summary',input$file1,Sys.Date(), '.txt')
    },
    content = function(file){
      # Load OMIM data
      omim_gene <- read.table("/home/anna/Desktop/sys_ndd_2025_subset.tsv", 
                        header=TRUE, sep="\t", stringsAsFactors=FALSE)

      # Join with statistical summary
      original_data <- stat_summ()
      joined_data <- merge(original_data, omim_gene, by="SYMBOL", all.x=TRUE)
      print(joined_data)

      # Write result
      write.table(joined_data, file, sep='\t', quote=F, row.names = F, col.names = T)
    }
  )

  # ============================================================================
  # OUTPUT: MYDATA TABLE (full dataset)
  # ============================================================================
  
  output$text <- DT::renderDataTable({
    validate(need(ncol(mydata()) != "0", "Please, upload your file"))
    mydata()
  })

  # ============================================================================
  # OUTPUT: FILTERED COVERAGE TABLE
  # ============================================================================
  
  output$text_cv <- DT::renderDataTable({
    validate(need(input$Gene_name != "" & input$Sample !="",
                  "Please select all required input: Gene, Chromosome,
                  Coverage threshold and Sample"))
    filtered_low()
  })

  # ============================================================================
  # OUTPUT: ALL GENE COVERAGE PLOT
  # ============================================================================
  
  output$all_gene <- renderPlot({
    validate(
      need(ncol(mydata()) != "0", "Unrecognized data set: Please
           upload your file"))
    options(shiny.sanitize.errors = TRUE)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Please wait a few minutes: Making plot",
                 detail = 'This may take a while', value = 0)
    for (i in 1:40) {
      progress$set(message = "Please wait a few minutes: Making plot",
                   detail = 'This may take a while', value = i)
      Sys.sleep(1.0)
    }
    Sys.sleep(0.1)
    p1()
  })

  # ============================================================================
  # OUTPUT: LOW COVERAGE POSITIONS TABLE
  # ============================================================================
  
  output$df.l <- DT::renderDataTable({
    table1()
  })

  # ============================================================================
  # OUTPUT: SEQUENCE PLOT
  # ============================================================================
  
  output$sequence <- renderPlot({
    validate(
      need(ncol(mydata()) != "0",
           "Unrecognized data set: Please load your file"))
    validate(
      need(input$Start_genomicPosition < input$end_genomicPosition,
           "Please select the right genomic position: end position is
           lower than start "))
    options(shiny.sanitize.errors = TRUE)
    progress <- shiny::Progress$new()
    for (i in 1:40) {
      progress$set(message = "Please wait a few minutes: Making plot",
                   detail = 'This may take a while', value = i)
      Sys.sleep(1.0)
    }
    on.exit(progress$close())
    Sys.sleep(0.1)
    p3()
  })

  # ============================================================================
  # OUTPUT: ANNOTATION SUMMARY TABLE
  # Conta varianti ClinVar e MutationAssessor H/M
  # ============================================================================
  
  output$tabExon <- DT::renderDataTable({
    validate(
      need(input$Gene_name != "" & input$Sample != "",
           "Please select all required input: Gene, Chromosome,
           Coverage threshold and Sample")
    )
    validate(
      need(try(!is.null(annotated_variants_data())),'Unrecognized coordinates:
      Please change exon number input and be sure that input box
      "Region coordinates" is filled.
           Please check your R environment has annotation file loaded.'))
    
    # USA annotated_variants_data() invece di intBED()
    data <- annotated_variants_data()
    
    if (is.null(data) || nrow(data) == 0) {
      return(NULL)
    }
    
    # Conta varianti importanti
    df <- data.frame(
      clinvar_count = length(which(data$ClinVar != '.')),
      high_medium_count = length(which(grepl("H|M", data$MutationAssessor)))
    )
    
    colnames(df) <- c('Low coverage positions in ClinVar',
                      'Low coverage positions with High or Medium impact')
    
    return(df)
  })

  # ============================================================================
  # OUTPUT: ANNOTATED VARIANTS TABLE (DT WIDGET)
  # TABELLA PRINCIPALE CON FORMATTING
  # ============================================================================
  
  output$uncover_position <- DT::renderDT({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Loading annotation table...", value = 0)
    Sys.sleep(0.1)
    
    # Validazione base
    validate(
      need(input$Gene_name != "" & input$Sample != "",
           "Please select Gene, Chromosome, Coverage threshold and Sample")
    )
    
    # Verifica che ci siano dati annotati
    data <- annotated_variants_data()
    
    if(is.null(data) || nrow(data) == 0) {
      showNotification("No annotated variants found. Check your annotation file path.", 
                       type = "warning", duration = 5)
      return(NULL)
    }
    
    cat("Displaying", nrow(data), "annotated variants\n")
    
    # Restituisci il widget DT
    condform_table()
  })

  # ============================================================================
  # DOWNLOAD: ANNOTATED VARIANTS EXCEL (con conditional formatting)
  # ============================================================================
  
  output$downloadData <- downloadHandler(
    filename = function() {
      sample_name <- if (!is.null(input$Sample)) {
        tools::file_path_sans_ext(input$Sample)
      } else {
        "variants"
      }
      paste0('annotated_variants_', sample_name, '_', Sys.Date(), '.xlsx')
    },
    
    content = function(file){
      # USA annotated_variants_data() invece di condform_table()!
      data <- annotated_variants_data()
      
      if (is.null(data) || nrow(data) == 0) {
        showNotification("No data to download!", type = "error")
        return(NULL)
      }
      
      # Rimuovi colonna helper prima di esportare
      data_export <- data
      data_export$highlight_important <- NULL
      
      # Crea workbook
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Variants")
      
      # Stili
      negStyle <- openxlsx::createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
      posStyle <- openxlsx::createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
      highlighted <- openxlsx::createStyle(fgFill = "yellow")
      hs <- openxlsx::createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF",
                        fontSize=12,
                        fontName="Arial Narrow", fgFill = "#4F80BD")
      
      # Scrivi dati
      openxlsx::writeData(wb, "Variants", data_export, headerStyle = hs)
      
      # Trova indici delle colonne dinamicamente
      col_names <- colnames(data_export)
      col_MutationAssessor <- which(col_names == "MutationAssessor")
      col_M_CAP <- which(col_names == "M_CAP")
      col_CADD_PHED <- which(col_names == "CADD_PHED")
      col_AF_gnomAD <- which(col_names == "AF_gnomAD")
      col_ClinVar <- which(col_names == "ClinVar")
      col_start <- which(col_names == "start")
      col_end <- which(col_names == "end")
      
      nrows <- nrow(data_export)
      
      # Conditional formatting: MutationAssessor
      if (length(col_MutationAssessor) > 0) {
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_MutationAssessor,
                              rows = 2:(nrows + 1), rule = '=="H"', style = negStyle)
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_MutationAssessor,
                              rows = 2:(nrows + 1), rule = '!="H"', style = posStyle)
      }
      
      # Conditional formatting: M_CAP
      if (length(col_M_CAP) > 0) {
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_M_CAP,
                              rows = 2:(nrows + 1), rule = '=="D"', style = negStyle)
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_M_CAP,
                              rows = 2:(nrows + 1), rule = '!="D"', style = posStyle)
      }
      
      # Conditional formatting: CADD_PHED
      if (length(col_CADD_PHED) > 0) {
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_CADD_PHED,
                              rows = 2:(nrows + 1), rule = ">=20", style = negStyle)
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_CADD_PHED,
                              rows = 2:(nrows + 1), rule = "<20", style = posStyle)
      }
      
      # Conditional formatting: AF_gnomAD
      if (length(col_AF_gnomAD) > 0) {
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_AF_gnomAD,
                              rows = 2:(nrows + 1), rule = "<=0.5", style = negStyle)
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_AF_gnomAD,
                              rows = 2:(nrows + 1), rule = ">0.5", style = posStyle)
      }
      
      # Conditional formatting: ClinVar
      if (length(col_ClinVar) > 0) {
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_ClinVar,
                              rows = 2:(nrows + 1), rule = '!="."', style = negStyle)
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_ClinVar,
                              rows = 2:(nrows + 1), rule = '=="."', style = posStyle)
      }
      
      # Highlight important variants (yellow background on start/end)
      if (length(col_start) > 0 && length(col_end) > 0) {
        important_rows <- which(
          grepl("H|M", data_export$MutationAssessor) & 
          data_export$ClinVar != "." & 
          !is.na(data_export$AF_gnomAD) & 
          data_export$AF_gnomAD < 0.5
        )
        
        if (length(important_rows) > 0) {
          openxlsx::addStyle(wb, "Variants", style = highlighted,
                             rows = important_rows + 1, 
                             cols = c(col_start, col_end),
                             gridExpand = TRUE)
        }
      }
      
      # Auto-size columns
      openxlsx::setColWidths(wb, "Variants", cols = 1:ncol(data_export), widths = "auto")
      
      # Salva
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      
      cat("Excel file saved:", nrow(data_export), "variants\n")
    }
  )

  # ============================================================================
  # OUTPUT: maxAF VALUE
  # ============================================================================
  
  output$maxAF <- renderText({
    maxAF_result <- data()
  
    if (is.null(maxAF_result) || is.na(maxAF_result[[1]])) {
      return("N/A - Check input parameters")
    }
  
    signif(maxAF_result[[1]], 3)
  })
  
  # ============================================================================
  # OUTPUT: maxAF TABLE (DT WIDGET con filtro maxAF)
  # ============================================================================
  
  output$uncoverPosition <- DT::renderDT({
    validate(
      need(ncol(mydata()) != "0",
           "Unrecognized data set: Please load your file"),
      need(input$Gene_name != "" &
             input$Sample != "", "Please select all required input: Gene,
           Chromosome, Coverage threshold and Sample")
    )
    
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Filtering by maxAF...", value = 0)
    Sys.sleep(0.1)

    # Verifica che ci sia almeno annotated_variants_data
    base_data <- annotated_variants_data()
    if(is.null(base_data) || nrow(base_data) == 0) {
      showNotification("No annotated data available for maxAF filtering!", 
                       type = "warning")
      return(NULL)
    }

    # Restituisci il widget DT (gestisce già il caso di 0 varianti con messaggio)
    uncover_maxaf()
  })

  # ============================================================================
  # DOWNLOAD: maxAF EXCEL
  # ============================================================================
  
  output$download_maxAF <- downloadHandler(
    filename = function() {
      paste0('uncover_maxAF_', Sys.Date(), '.xlsx')
    },
    content = function(file){
      # USA il dataframe pulito, NON il widget!
      data <- uncover_maxaf_data()
      
      if (is.null(data) || nrow(data) == 0) {
        showNotification("No data to download!", type = "error")
        return(NULL)
      }
      
      # Rimuovi colonne helper
      data_export <- data
      data_export$highlight_maxaf <- NULL
      data_export$highlight_important <- NULL
      
      wb1 <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb1, "MaxAF")
      
      negStyle <- openxlsx::createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
      posStyle <- openxlsx::createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
      highlighted <- openxlsx::createStyle(fgFill = "yellow")
      hs <- openxlsx::createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF",
                        fontSize=12,
                        fontName="Arial Narrow", fgFill = "#4F80BD")
      
      openxlsx::writeData(wb1, "MaxAF", data_export, headerStyle = hs)
      
      # Trova indici colonne dinamicamente
      col_names <- colnames(data_export)
      col_MutationAssessor <- which(col_names == "MutationAssessor")
      col_M_CAP <- which(col_names == "M_CAP")
      col_CADD_PHED <- which(col_names == "CADD_PHED")
      col_AF_gnomAD <- which(col_names == "AF_gnomAD")
      col_ClinVar <- which(col_names == "ClinVar")
      col_start <- which(col_names == "start")
      col_end <- which(col_names == "end")
      
      nrows <- nrow(data_export)
      
      # Conditional formatting (stesso pattern di downloadData)
      if (length(col_MutationAssessor) > 0) {
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_MutationAssessor,
                              rows = 2:(nrows + 1), rule = '=="H"', style = negStyle)
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_MutationAssessor,
                              rows = 2:(nrows + 1), rule = '!="H"', style = posStyle)
      }
      
      if (length(col_M_CAP) > 0) {
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_M_CAP,
                              rows = 2:(nrows + 1), rule = '=="D"', style = negStyle)
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_M_CAP,
                              rows = 2:(nrows + 1), rule = '!="D"', style = posStyle)
      }
      
      if (length(col_CADD_PHED) > 0) {
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_CADD_PHED,
                              rows = 2:(nrows + 1), rule = ">=20", style = negStyle)
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_CADD_PHED,
                              rows = 2:(nrows + 1), rule = "<20", style = posStyle)
      }
      
      if (length(col_AF_gnomAD) > 0) {
        # Per maxAF usiamo threshold 0.05 invece di 0.5
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_AF_gnomAD,
                              rows = 2:(nrows + 1), rule = "<=0.05", style = negStyle)
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_AF_gnomAD,
                              rows = 2:(nrows + 1), rule = ">0.05", style = posStyle)
      }
      
      if (length(col_ClinVar) > 0) {
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_ClinVar,
                              rows = 2:(nrows + 1), rule = '!="."', style = negStyle)
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_ClinVar,
                              rows = 2:(nrows + 1), rule = '=="."', style = posStyle)
      }
      
      # Highlight important variants
      if (length(col_start) > 0 && length(col_end) > 0) {
        maxAF_value <- signif(data()[[1]], 3)
        
        important_rows <- which(
          grepl("H|M", data_export$MutationAssessor) & 
          data_export$ClinVar != "." & 
          !is.na(data_export$AF_gnomAD) & 
          data_export$AF_gnomAD < maxAF_value
        )
        
        if (length(important_rows) > 0) {
          openxlsx::addStyle(wb1, "MaxAF", style = highlighted,
                             rows = important_rows + 1, 
                             cols = c(col_start, col_end),
                             gridExpand = TRUE)
        }
      }
      
      openxlsx::setColWidths(wb1, "MaxAF", cols = 1:ncol(data_export), widths = "auto")
      openxlsx::saveWorkbook(wb1, file, overwrite = TRUE)
      
      cat("MaxAF Excel saved:", nrow(data_export), "variants\n")
    }
  )
  
  
  # ============================================================================
  # CLOSE APP BUTTON
  # ============================================================================
  
  observe({
    if (input$close > 0) stopApp()
  })
}