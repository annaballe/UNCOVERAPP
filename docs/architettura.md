# Architettura e struttura del codice

Struttura principale del pacchetto:
- R/: funzioni esportate (run.uncoverapp, buildInput, annotate_variants, annotate_all_lowcov, getAnnotationFiles)
- inst/shiny-dir/: app Shiny (ui.R, server.R, compute-*.R, assets in www/)
- inst/extdata/: markdown introduttivi e file di esempio
- vignettes/: vignetta R Markdown
- man/: documentazione Rd generata da roxygen2
- tests/: testthat (unit) e inst/shiny-dir/tests (test di app)
- DESCRIPTION, NAMESPACE: metadati pacchetto

Shiny app (inst/shiny-dir):
- ui.R: definizione interfaccia (tab Home, Processing and Statistical Summary, Annotation, Plots, Tables, Binomial, Max AF)
- server.R: logica server e orchestrazione
- compute-*.R: moduli di calcolo (preprocess, annotation, plots, tables, binomial, maxAF)

Funzioni chiave R:
- run.uncoverapp/ uncoverAPP: risoluzione del percorso app e lancio con shiny::runApp
- buildInput: pipeline per generare BED di copertura (da BAM o BED) e sommari statistici
- annotate_variants: annotazione low coverage in un gene/regione
- annotate_all_lowcov: annotazione low coverage su tutto il dataset
- getAnnotationFiles: gestione cache BiocFileCache per file dbNSFP
