suppressPackageStartupMessages({
        require(OrganismDbi)
        require(TxDb.Hsapiens.UCSC.hg19.knownGene)
        require(TxDb.Hsapiens.UCSC.hg38.knownGene)
        require(org.Hs.eg.db)
        require(rlist)
        require(Rsamtools)
    })
    cat("Heavy packages loaded!\n")



ff <- reactiveValues()

gene_list <- reactive({
   file_gene <- input$gene1
   if (is.null(input$gene1)) return(NULL)
   
   if(input$type_file == "target_bed") {
      tmp_gene1 <- read.table(input$gene1$datapath, stringsAsFactors = F)
      tmp_gene <- tmp_gene1[1:4]
      colnames(tmp_gene) <- c('chr', 'start', 'end', 'SYMBOL')
      ff <- tmp_gene
   } else {
      tmp_gene <- scan(file_gene$datapath, character(), quote = "")
      ff <- tmp_gene
   }
})

list_coverage <- reactive({
    if (is.null(input$bam_list)) return(NULL)
    scan(input$bam_list$datapath, character(), quote = "")
})

all_gene <- reactive({
   if (input$Genome == "hg19") {
      TxDb.Hsapiens.UCSC.hg19.knownGene
   } else {
      TxDb.Hsapiens.UCSC.hg38.knownGene
   }
})

# ============================================================================
# NUOVO REACTIVE CENTRALIZZATO: calcola UNA SOLA VOLTA la validazione dei geni
# ============================================================================
validated_genes_data <- reactive({
   if (is.null(gene_list())) return(NULL)
   if (is.null(list_coverage())) return(NULL)
   if (input$type_file == "target_bed") return(NULL)  # Non serve per bed files
   
   cat("\n=== GENE VALIDATION (CENTRALIZED) ===\n")
   
   # STAGE 1: Convert gene names to ENTREZ IDs (UNA SOLA VOLTA!)
   cat("Stage 1: Converting gene names to ENTREZ IDs...\n")
   
   my_gene_name <- OrganismDbi::select(
      org.Hs.eg.db, 
      key = gene_list(), 
      columns = "ENTREZID", 
      keytype = "ALIAS"
   ) %>%
      tibble::as_tibble() %>%
      dplyr::filter(!is.na(ENTREZID)) %>%
      dplyr::group_by(ALIAS) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
   
   found_genes_stage1 <- unique(my_gene_name$ALIAS)
   not_found_stage1 <- base::setdiff(gene_list(), found_genes_stage1)
   
   if (length(not_found_stage1) > 0) {
      cat("WARNING: Genes not found in org.Hs.eg.db:\n")
      for(gene in not_found_stage1) {
         cat(paste("  - Gene '", gene, "'\n", sep=""))
      }
      cat("\n")
   }
   
   if (nrow(my_gene_name) == 0) {
      stop("ERROR: No valid genes found in org.Hs.eg.db.")
   }
   
   cat(paste("Stage 1 result:", nrow(my_gene_name), "mappings from", 
             length(found_genes_stage1), "genes\n\n"))
   
   # STAGE 2: Validate against genome (UNA SOLA VOLTA!)
   cat("Stage 2: Validating against reference genome...\n")
   
   all_entrez_ids <- unique(my_gene_name$ENTREZID)
   genome_keys <- AnnotationDbi::keys(all_gene(), "GENEID")  # CACHE THIS!
   valid_entrez_ids <- base::intersect(all_entrez_ids, genome_keys)
   invalid_entrez_ids <- base::setdiff(all_entrez_ids, genome_keys)
   
   invalid_genes_stage2 <- character(0)
   if (length(invalid_entrez_ids) > 0) {
      invalid_genes_stage2 <- my_gene_name %>% 
         dplyr::filter(ENTREZID %in% invalid_entrez_ids) %>%
         dplyr::pull(ALIAS) %>%
         unique()
      
      cat("WARNING: Genes not in reference genome:\n")
      for(gene in invalid_genes_stage2) {
         cat(paste("  - Gene '", gene, "' (", input$Genome, ")\n", sep=""))
      }
      cat("\n")
   }
   
   # Filter to valid genes only
   my_gene_name <- my_gene_name %>%
      dplyr::filter(ENTREZID %in% valid_entrez_ids)
   
   if (length(valid_entrez_ids) == 0) {
      stop("ERROR: No genes passed validation.")
   }
   
   cat(paste("Stage 2 result:", length(valid_entrez_ids), "validated genes\n"))
   
   total_excluded <- length(not_found_stage1) + length(invalid_genes_stage2)
   if (total_excluded > 0) {
      cat(paste("SUMMARY: ", total_excluded, " genes excluded (", 
                length(not_found_stage1), " stage 1, ", 
                length(invalid_genes_stage2), " stage 2)\n", sep=""))
   }
   cat("=== VALIDATION COMPLETE ===\n\n")
   
   # RETURN: tutto quello che serve ai reactive successivi
   list(
      my_gene_name = my_gene_name,
      valid_entrez_ids = valid_entrez_ids,
      genome_keys = genome_keys,  # Cache per riuso
      not_found_stage1 = not_found_stage1,
      invalid_genes_stage2 = invalid_genes_stage2
   )
})

# ============================================================================
# SEMPLIFICATO: usa solo i risultati pre-calcolati
# ============================================================================
no_entrID <- reactive({
   if (is.null(gene_list())) return(NULL)
   if (is.null(list_coverage())) return(NULL)
   
   if (input$type_file == "target_bed") {
      return(data.frame(matrix(ncol = 2, nrow = 0)))
   }
   
   # Usa i dati già validati!
   validated <- validated_genes_data()
   
   if (is.null(validated)) {
      return(data.frame(matrix(ncol = 2, nrow = 0)))
   }
   
   # Se ci sono geni validi, ritorna dataframe vuoto (nessun errore)
   if (length(validated$valid_entrez_ids) > 0) {
      return(data.frame(matrix(ncol = 2, nrow = 0)))
   }
   
   # Altrimenti ritorna errori
   errore <- data.frame(
      Gene = c("No valid genes found"),
      Issue = c("All genes failed validation"),
      stringsAsFactors = FALSE
   )
   return(errore)
})

# ============================================================================
# OTTIMIZZATO: usa dati pre-calcolati + ottimizzazioni loop
# ============================================================================
for_bed <- reactive({
   if (is.null(gene_list())) return(NULL)
   if (is.null(list_coverage())) return(NULL)
   if (nrow(no_entrID()) != 0) return(no_entrID())
   
   if (input$type_file == "target_bed") {
      return(gene_list())
   }
   
   # USA I DATI GIÀ VALIDATI - zero duplicazione!
   validated <- validated_genes_data()
   if (is.null(validated)) return(NULL)
   
   my_gene_name <- validated$my_gene_name
   ID <- validated$valid_entrez_ids
   # genome_keys già disponibile ma non serve più (vedi sotto)
   
   # OTTIMIZZAZIONE: carica exonsBy UNA SOLA VOLTA (era dentro il loop!)
   cds <- OrganismDbi::exonsBy(all_gene(), by = "tx", use.names = TRUE)
   
   # OTTIMIZZAZIONE: usa lapply invece di loop con rbind
   pre_list <- lapply(ID, function(i) {
      # Rimuovo il check ridondante: if (i %in% genome_keys)
      # È già garantito che i è valido!
      
      txid <- OrganismDbi::select(all_gene(), i, "TXNAME", "GENEID")[["TXNAME"]]
      
      if (length(txid) == 0) return(NULL)
      
      coor <- as.data.frame(cds[names(cds) %in% txid])
      
      if (nrow(coor) == 0) return(NULL)
      
      coor$ENTREZID <- i
      return(coor)
   })
   
   # Combina tutto in una volta (molto più veloce di rbind iterativo)
   pre <- dplyr::bind_rows(pre_list) %>%
      dplyr::inner_join(my_gene_name, by = "ENTREZID") %>%
      dplyr::mutate(
         start = start - 10,
         end = end + 10
      )
   
   # OTTIMIZZAZIONE: elimina il secondo loop, usa operazioni vettoriali
   for_bed <- pre %>%
      dplyr::transmute(
         chr = as.character(seqnames),
         start = start,
         end = end,
         SYMBOL = ALIAS
      ) %>%
      dplyr::distinct()
   
   # Apply notation transformation
   if (input$notation == "number") {
      for_bed <- for_bed %>%
         dplyr::mutate(chr = sub("^.{0,3}", "", chr, perl = TRUE))
   }
   
   return(for_bed)
})

# ============================================================================
# Il resto del codice rimane invariato
# ============================================================================
coverage_input <- reactive({
    if (is.null(gene_list())) return(NULL)
    if (!is.null(no_entrID()) && nrow(no_entrID()) != 0) return(no_entrID())
    if (is.null(list_coverage())) return(NULL)

    if (input$type_coverage == "bam") {
        for_grange <- GenomicRanges::makeGRangesFromDataFrame(for_bed(), keep.extra.columns = TRUE)
        param <- Rsamtools::ScanBamParam(which=for_grange)

        p_param <- Rsamtools::PileupParam(
            distinguish_nucleotides=FALSE,
            distinguish_strands=FALSE,
            min_mapq=as.numeric(input$MAPQ),
            min_base_quality=as.numeric(input$base_qual),
            min_nucleotide_depth=1,
            max_depth=150000
        )

        df <- list()
        for (i in list_coverage()) {
            pileup_df <- Rsamtools::pileup(i, scanBamParam=param, pileupParam=p_param)
            df <- rlist::list.append(df, pileup_df)
        }

        lst1 <- lapply(df, function(x) transform(x[,-5]))
        lst2 <- lapply(lst1, function(x) transform(x[!duplicated(x),]))
        riarrange.df <- function(list_df) {
            list_df %>%
                dplyr::mutate(end=pos) %>%
                dplyr::group_by(seqnames,pos,end) %>%
                dplyr::summarise(count=sum(count)) %>%
                dplyr::arrange()
        }
        lst3 <- lapply(lst2, riarrange.df)
        pp<- Reduce(function(...) merge(..., by=c("seqnames","pos","end")), lst3)
        pp[is.na(pp)] <- 0
        colnames(pp)[1:3] <- c("seqnames","start","end")
        pp$seqnames <- paste0("chr", sub("^chr", "", as.character(pp$seqnames)))
        return(pp)
    }
    
    if (input$type_coverage == "bed") {
        for_grange <- GenomicRanges::makeGRangesFromDataFrame(
            for_bed(), 
            keep.extra.columns = TRUE
        )
        
        df_list <- list()
        for (file in list_coverage()) {
            df <- read.table(file, 
                            header = FALSE,
                            stringsAsFactors = FALSE,
                            colClasses = c("character", "integer", "integer", "numeric"))
            
            colnames(df) <- c("seqnames", "start", "end", "count")
            df_gr <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
            df_gr_filtered <- IRanges::subsetByOverlaps(df_gr, for_grange, type="any")
            df_filtered <- as.data.frame(df_gr_filtered)
            df_filtered <- df_filtered[, c("seqnames", "start", "end", "count")]
            df_list <- append(df_list, list(df_filtered))
        }
        
        pp <- Reduce(function(...) merge(..., by=c("seqnames","start","end"), all=TRUE), df_list)
        pp[is.na(pp)] <- 0
        pp$seqnames <- factor(pp$seqnames)
        pp$start <- as.integer(pp$start)
        pp$end <- as.integer(pp$end)
        count_cols <- grep("^count", colnames(pp), value = TRUE)
        pp[count_cols] <- lapply(pp[count_cols], as.integer)
        pp$seqnames <- paste0("chr", sub("^chr", "", as.character(pp$seqnames)))
        
        return(pp)
    }
})

name_sample <- reactive({
   if (is.null(list_coverage())) return(character(0))
   tools::file_path_sans_ext(basename(list_coverage()))
})

stat_summ <- reactive({
   if (is.null(gene_list())) return(NULL)
   if (is.null(list_coverage())) return(NULL)
   if (nrow(no_entrID()) != 0) return(no_entrID())
   
   ppinp <- as.data.frame(coverage_input())
   colnames(ppinp)[1:3] <- c("seqnames","start","end")
   n <- length(colnames(ppinp)[-1:-3])
   samples <- name_sample()
   colnames(ppinp)[-1:-3] <- paste0("sample_", samples)[seq_len(n)]

   for_grange <- GenomicRanges::makeGRangesFromDataFrame(for_bed(),
                                                         keep.extra.columns = TRUE)

   for_range_pp <- GenomicRanges::makeGRangesFromDataFrame(ppinp,
                                                           keep.extra.columns = TRUE)
   tp <- GenomicRanges::findOverlaps(query= for_range_pp,
                                    subject = for_grange, type="any",
                                    select = "all")
   sts_df <- data.frame(for_range_pp[queryHits(tp),], for_grange[subjectHits(tp),])

   statistiche <- sts_df[!duplicated(sts_df$start),]
   statistiche <- subset(statistiche,
                        select = -c(width, strand, seqnames.1,
                                    start.1, end.1, width.1, strand.1))

   colnames(statistiche)[1:3] <- c("chromosome","start","end")
   merge_g <- dplyr::full_join(for_bed(),statistiche, by="SYMBOL")
   
   col_name <- colnames(merge_g)
   col.sub <- col_name[grepl("sample_", col_name)]
   merge_g[col.sub] <- sapply(merge_g[col.sub],as.numeric)

   x <- list()
   for (i in col.sub){
      x_df <- merge_g %>%
         dplyr::select("chromosome","start.x","end.x", "SYMBOL",i)
      colnames(x_df)[5] <- "value"
      x_df$sample <- paste0(i)
      x <- list.append(x, x_df)
   }
   
   statistical_operation <- function(df){
      df %>%
         dplyr::group_by(SYMBOL,sample) %>%
         dplyr::summarize(Total= sum(value, na.rm=TRUE),
                          Mean = mean(value, na.rm=TRUE),
                          Median= median(value, na.rm=TRUE),
                          number_of_position_under_20x = sum(value < 20),
                          percentage_under_20x= (sum(value < 20)/sum(value, na.rm=TRUE))*100) %>%
         as.data.frame(.) %>%
         dplyr::mutate_if(is.numeric, round, 3)
   }

   out_r <- do.call(rbind, lapply(x, function(x) statistical_operation(x)))
   return(out_r)
})