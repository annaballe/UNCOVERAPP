# Installazione

Prerequisiti:
- R >= 4.0.0
- Java installato (per alcuni pacchetti a corredo)
- File di annotazione dbNSFP preparati (sorted.bed.gz e sorted.bed.gz.tbi)

Installazione da GitHub (sviluppo):
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
install.packages("devtools")
devtools::install_github("Manuelaio/uncoverappLib")
```

Scarico file di annotazione in cache locale (Zenodo):
```r
library(uncoverappLib)
getAnnotationFiles(verbose = TRUE)
```

Variabili d’ambiente opzionali per scegliere file di annotazione custom:
- UNCOVERAPP_HG19_ANNOTATION=/percorso/sorted_hg19.bed.gz
- UNCOVERAPP_HG38_ANNOTATION=/percorso/sorted_hg38.bed.gz
