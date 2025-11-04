# API delle funzioni R

## run.uncoverapp(where = c("browser","viewer","window"))
- Avvia l’app Shiny in browser/Viewer/nuova finestra RStudio.

## uncoverAPP()
- Risolve il percorso dell’app in inst/shiny-dir e lancia shiny::runApp.

## buildInput(geneList, genome, type_bam, bamList, outDir, type_input,
             MAPQ.min=1, base.quality=1, type_coverage="bam",
             input_coord_system="1-based", annotation_file=NULL)
- Input: lista geni o target BED; lista file coverage (BAM/BED); parametri genome/coordinate.
- Output: nella cartella outDir/output/ due file
  - <data>.bed: copertura per intervallo/posizione (colonne: chromosome,start,end,count_<campioni>)
  - <data>_statistical_summary.txt: sommario per gene x campione (media, mediana, %<20x, ecc.)

## annotate_variants(sample_data, target_sample, gene,
                     coverage_threshold=20, query_region="", genome="hg19",
                     annotation_file=NULL, output_intersect, output_formatted)
- Filtra posizioni <= soglia nel gene/regione; interseca con dbNSFP via Tabix; esporta TSV ed Excel con formattazioni.
- Gestione flessibile del nome colonna campione: "<nome>", "count_<nome>", match parziale.

## annotate_all_lowcov(sample_data, target_sample,
                       coverage_threshold=20, genome="hg19",
                       annotation_file=NULL, output_intersect, output_formatted)
- Come sopra, ma genome-wide (non limitato a un gene/regione).

## getAnnotationFiles(verbose=FALSE, vignette=FALSE)
- Scarica in cache BiocFileCache i file sorted.bed.gz e .tbi da Zenodo e restituisce i percorsi.
