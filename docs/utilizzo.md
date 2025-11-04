# Utilizzo

Avvio dell’app Shiny (GUI):
```r
library(uncoverappLib)
run.uncoverapp(where = "window")  # "browser" | "viewer" | "window"
```

Formati di input accettati dalla GUI:
- Lista geni (.txt): un nome HGNC per riga
- Lista coverage (.list): percorsi assoluti a file coverage (BAM o BED), uno per riga

Workflow da riga di comando (senza GUI):
1) Preparazione input di copertura:
```r
buildInput(
  geneList = "genes.txt" | "target.bed",
  genome = "hg19" | "hg38",
  type_bam = "chr" | "number",
  bamList = "coverage.list",
  outDir = "./output",
  type_input = "genes" | "target",
  type_coverage = "bam" | "bed",
  input_coord_system = "1-based" | "0-based",
  annotation_file = NULL
)
```
2) Annotazione per gene o regione specifica:
```r
annotate_variants(
  sample_data = "./output/<data>.bed",
  target_sample = "<nome_campione>" | "count_<nome>",
  gene = "BRCA1",
  coverage_threshold = 20,
  query_region = "chr10:73697750-73697765",
  genome = "hg19",
  annotation_file = "/path/sorted_hg19.bed.gz"
)
```
3) Annotazione di tutte le posizioni low coverage (genome-wide):
```r
annotate_all_lowcov(
  sample_data = "./output/<data>.bed",
  target_sample = "<nome_campione>",
  coverage_threshold = 20,
  genome = "hg19",
  annotation_file = "/path/sorted_hg19.bed.gz"
)
```
