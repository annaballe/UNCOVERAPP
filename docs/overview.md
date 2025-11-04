# Panoramica

unCOVERApp è una web-app Shiny per:
- Evidenziare posizioni genomiche con copertura insufficiente rispetto a una soglia definita dall’utente
- Annotare tali posizioni con dbNSFP (prediction scores, AF gnomAD, ClinVar/OMIM, ecc.)
- Fornire strumenti statistici di supporto (probabilità binomiale di adeguata copertura; Maximum Credible Allele Frequency)

Moduli principali dell’app:
- Processing and Statistical Summary: preparazione input e sommari per gene/campione
- Annotation: annotazione dbNSFP delle posizioni low coverage
- Plots e Tables: visualizzazioni interattive
- Binomial: calcolo probabilità di copertura
- Max AF: calcolo Maximum Credible Allele Frequency

Flusso tipico:
1) Scarico dei file di annotazione (getAnnotationFiles)
2) Preparazione input (buildInput) a partire da lista geni o target BED e lista coverage (BAM/BED)
3) Avvio GUI (run.uncoverapp) e analisi/visualizzazione
4) Opzionale: annotazione batch via funzioni R (annotate_variants, annotate_all_lowcov)
