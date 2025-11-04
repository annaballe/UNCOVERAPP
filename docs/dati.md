# Dati e file di annotazione

Annotazione dbNSFP:
- sorted.bed.gz: BED bgzip/tabix con colonne selezionate da dbNSFP v4.0
- sorted.bed.gz.tbi: indice Tabix
- Download automatico via getAnnotationFiles() su Zenodo

Esempi e risorse (inst/extdata/):
- intro.md, prep_input.md, CONTACTS.md
- Esempi: POLG.example.bed, example_POLG.bam (.bai), mygene.txt

Input supportati:
- Geni: .txt (un gene per riga)
- Target: BED con colonne chr, start, end, SYMBOL
- Coverage: BAM (pileup) o BED di coverage 0/1-based (specificare input_coord_system)

Output principali:
- coverage .bed con colonne chromosome,start,end,count_<sample>
- statistiche per gene/campione (.txt)
- annotazioni: TSV ed Excel formattato
