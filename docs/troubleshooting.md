# Risoluzione problemi

- File di annotazione non trovati
  - Usare getAnnotationFiles() oppure impostare UNCOVERAPP_HG19_ANNOTATION / UNCOVERAPP_HG38_ANNOTATION.
  - Verificare anche la presenza dell’indice .tbi.

- Colonna campione non riconosciuta in annotate_*:
  - Provare sia "<nome>" che "count_<nome>"; è supportato anche il match parziale case-insensitive.

- Geni non riconosciuti in buildInput:
  - Controllare log in outDir/output/ (preprocessing_log*.txt). Usare simboli HGNC o fornire un BED di target.

- Disallineamento notazione cromosomi (chr vs numero):
  - Impostare type_bam coerente ("chr" o "number").

- Coordinate 0-based vs 1-based nei BED:
  - Specificare input_coord_system correttamente per evitare shift di 1 base.
