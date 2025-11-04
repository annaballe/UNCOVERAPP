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
   start_time <- Sys.time()
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
    cat(paste("Load input files:", round(difftime(Sys.time(), start_time, units="mins"), 2), "minutes\n\n"))
    return(ff)
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
   start_time <- Sys.time()
   cat("\n=== RUNNING ON R SESSION:", Sys.info()["nodename"], "===\n")
   cat("R Process ID:", Sys.getpid(), "\n")
   cat("Working directory:", getwd(), "\n\n")

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
   cat(paste("Gene validation took:", round(difftime(Sys.time(), start_time, units="mins"), 2), "minutes\n\n"))
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
   start_time <- Sys.time()
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
   txid_map <- OrganismDbi::select(all_gene(), ID, "TXNAME", "GENEID")
   
   # OTTIMIZZAZIONE: usa lapply invece di loop con rbind
    pre_list<- lapply(split(txid_map, txid_map$GENEID), function(gene_txs) {
        tx_names <- gene_txs$TXNAME
        coor <- as.data.frame(cds[names(cds) %in% tx_names])
        if (nrow(coor) == 0) return(NULL)
        coor$ENTREZID <- gene_txs$GENEID[1]
        return(coor)
        })

        pre_list <- pre_list[!sapply(pre_list, is.null)]
    
    cat("DEBUG: pre_list has", length(pre_list), "non-NULL elements\n")
    
    # Check se ci sono dati
    if (length(pre_list) == 0) {
        stop("No valid transcript coordinates found for any gene!")
    }
   
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
   cat(paste("BED creation took:", round(difftime(Sys.time(), start_time, units="mins"), 2), "minutes\n\n"))

   
   return(for_bed)
})

# ============================================================================
# COMPLETELY FIXED COVERAGE PROCESSING
# Handles both 0-based and 1-based coordinate systems
# Coverage BED format: chr start end depth (uniform depth per interval)
# ============================================================================

coverage_input <- reactive({
    if (is.null(gene_list())) return(NULL)
    if (!is.null(no_entrID()) && nrow(no_entrID()) != 0) return(no_entrID())
    if (is.null(list_coverage())) return(NULL)
    
    coverage_time <- Sys.time()
    cat("\n=== COVERAGE PROCESSING ===\n")
    
    # Determine input coordinate system (add this as UI input!)
    # For now, default to 0-based (standard BED format)
    input_coord_system <- if (!is.null(input$input_coordinate_system)) {
        input$input_coordinate_system
    } else {
        "0-based"  # Default assumption
    }
    cat(paste("Input coordinate system:", input_coord_system, "\n"))
    
    # Get target regions and convert to 1-based for GenomicRanges
    target_bed <- for_bed()
    
    for_grange <- GenomicRanges::makeGRangesFromDataFrame(
        target_bed, 
        keep.extra.columns = TRUE
    )
    
    if (input$type_coverage == "bam") {
        # BAM processing (pileup returns per-base 1-based positions)
        cat("Processing BAM files with pileup...\n")
        
        param <- Rsamtools::ScanBamParam(which = for_grange)
        p_param <- Rsamtools::PileupParam(
            distinguish_nucleotides = FALSE,
            distinguish_strands = FALSE,
            min_mapq = as.numeric(input$MAPQ),
            min_base_quality = as.numeric(input$base_qual),
            min_nucleotide_depth = 1,
            max_depth = 150000
        )
        
        lst3 <- lapply(list_coverage(), function(i) {
            pileup_df <- Rsamtools::pileup(i, scanBamParam = param, pileupParam = p_param)
            pileup_df <- pileup_df[, -5]  # Remove nucleotide column
            pileup_df <- pileup_df[!duplicated(pileup_df), ]
            
            # Pileup returns 1-based positions
            pileup_df %>%
                dplyr::mutate(end = pos) %>%
                dplyr::group_by(seqnames, pos, end) %>%
                dplyr::summarise(count = sum(count), .groups = "drop")
        })
        
        pp <- Reduce(function(...) merge(..., by = c("seqnames", "pos", "end")), lst3)
        pp[is.na(pp)] <- 0
        colnames(pp)[1:3] <- c("seqnames", "start", "end")
        
        cat("BAM pileup: per-base coverage (1-based positions)\n")
        
    } else if (input$type_coverage == "bed") {
        cat("Processing BED coverage files (format: chr start end depth)...\n")
        cat("Each interval represents a uniform coverage block\n\n")
        
        # Read and process all coverage BED files
        df_list <- lapply(seq_along(list_coverage()), function(idx) {
            file <- list_coverage()[idx]
            cat(paste("Sample", idx, ":", basename(file), "\n"))
            
            df <- read.table(
                file, 
                header = FALSE,
                stringsAsFactors = FALSE,
                colClasses = c("character", "integer", "integer", "numeric")
            )
            colnames(df) <- c("seqnames", "start", "end", "count")
            
            cat(paste("  Raw intervals:", nrow(df), "\n"))
            
            # Convert to 1-based if input is 0-based
            if (input_coord_system == "0-based") {
                df <- df %>% dplyr::mutate(start = start + 1)
                cat("  Converted from 0-based to 1-based\n")
            }
            
            # CRITICAL FIX: Add unique identifier to preserve count through intersect
            df$interval_id <- paste0("interval_", seq_len(nrow(df)))
            
            # Convert to GRanges (now definitely 1-based, with metadata preserved)
            df_gr <- GenomicRanges::makeGRangesFromDataFrame(
                df, 
                keep.extra.columns = TRUE  # Keeps count AND interval_id
            )
            
            # Find overlaps with target regions
            overlaps <- GenomicRanges::findOverlaps(df_gr, for_grange, type = "any")
            
            if (length(overlaps) == 0) {
                warning(paste("No overlaps found for", basename(file)))
                return(data.frame(
                    seqnames = character(), 
                    start = integer(), 
                    end = integer(), 
                    count = numeric(),
                    stringsAsFactors = FALSE
                ))
            }
            
            cat(paste("  Overlapping intervals:", length(unique(queryHits(overlaps))), "\n"))
            
            # CRITICAL FIX: Use pintersect to trim AND preserve metadata
            # pintersect trims each query interval to the overlapping subject intervals
            df_gr_overlapping <- df_gr[queryHits(overlaps)]
            for_grange_overlapping <- for_grange[subjectHits(overlaps)]
            
            # pintersect automatically preserves all metadata (count, interval_id)
            df_gr_trimmed <- GenomicRanges::pintersect(
                df_gr_overlapping,
                for_grange_overlapping,
                ignore.strand = TRUE
            )
            
            cat(paste("  Trimmed to", length(df_gr_trimmed), "intervals\n"))
            
            # Convert back to dataframe - count is automatically preserved!
            df_trimmed <- as.data.frame(df_gr_trimmed)
            
            # Clean up: remove metadata columns we don't need
            df_trimmed <- df_trimmed %>%
                dplyr::select(seqnames, start, end, count) %>%
                dplyr::distinct()  # Remove any duplicate intervals created by overlapping targets
            
            cat(paste("  Final intervals:", nrow(df_trimmed), "\n\n"))
            
            return(df_trimmed)
        })
        
        # Remove empty dataframes
        df_list <- df_list[sapply(df_list, nrow) > 0]
        
        if (length(df_list) == 0) {
            stop("ERROR: No coverage data overlaps with target regions!")
        }
        
        cat(paste("Merging", length(df_list), "trimmed coverage files...\n"))
        
        # Merge all samples by genomic coordinates
        if (length(df_list) == 1) {
            pp <- df_list[[1]]
        } else {
            pp <- df_list[[1]]
            for (i in 2:length(df_list)) {
                pp <- merge(
                    pp, 
                    df_list[[i]], 
                    by = c("seqnames", "start", "end"), 
                    all = TRUE,
                    suffixes = c("", paste0(".sample", i))
                )
            }
        }
        
        # Replace NA with 0 for missing coverage
        pp[is.na(pp)] <- 0
        
        cat(paste("Merged result:", nrow(pp), "unique intervals\n"))
    }
    
    # Normalize chromosome names
    pp$seqnames <- paste0("chr", sub("^chr", "", as.character(pp$seqnames)))
    
    cat(paste("\nFinal coverage dimensions:", nrow(pp), "intervals x", ncol(pp)-3, "samples\n"))
    cat(paste("Coverage processing took:", 
              round(difftime(Sys.time(), coverage_time, units = "mins"), 2), 
              "minutes\n"))
    
    # Store coordinate system metadata
    attr(pp, "coordinate_system") <- "1-based"
    attr(pp, "trimmed_to_targets") <- TRUE
    attr(pp, "original_input_system") <- input_coord_system
    
    return(pp)
})

name_sample <- reactive({
   if (is.null(list_coverage())) return(character(0))
   tools::file_path_sans_ext(basename(list_coverage()))
})

# ============================================================================
# COMPLETELY FIXED STATISTICS with proper weighted calculations
# ============================================================================

stat_summ <- reactive({
    start_time <- Sys.time()
    if (is.null(gene_list())) return(NULL)
    if (is.null(list_coverage())) return(NULL)
    if (!is.null(no_entrID()) && nrow(no_entrID()) != 0) return(no_entrID())

    cat("\n=== STATISTICS CALCULATION ===\n")
    
    # Use 1-based coverage for statistics
    ppinp <- as.data.frame(coverage_input())
    colnames(ppinp)[1:3] <- c("seqnames", "start", "end")
    
    # Identify sample columns
    n <- length(colnames(ppinp)[-1:-3])
    samples <- name_sample()
    
    # Ensure we have the right number of samples
    if (length(samples) != n) {
        warning(paste("Mismatch: found", n, "data columns but", length(samples), "sample names"))
        samples <- paste0("sample_", seq_len(n))
    }
    
    colnames(ppinp)[-1:-3] <- paste0("sample_", samples)[seq_len(n)]

    # Target regions (convert to 1-based if needed)
    target_bed <- for_bed()
    input_coord_system <- if (!is.null(input$input_coordinate_system)) {
        input$input_coordinate_system
    } else {
        "0-based"
    }
    
    if (input_coord_system == "0-based") {
        target_bed <- target_bed %>% dplyr::mutate(start = start + 1)
    }
    
    for_grange <- GenomicRanges::makeGRangesFromDataFrame(
        target_bed,
        keep.extra.columns = TRUE
    )

    for_range_pp <- GenomicRanges::makeGRangesFromDataFrame(
        ppinp,
        keep.extra.columns = TRUE
    )

    # Find overlaps between coverage and targets
    tp <- GenomicRanges::findOverlaps(
        query = for_range_pp,
        subject = for_grange,
        type = "any",
        select = "all"
    )

    # Merge coverage with gene annotations
    sts_df <- data.frame(
        for_range_pp[queryHits(tp), ],
        for_grange[subjectHits(tp), ]
    )

    # Remove duplicates and clean up columns
    statistiche <- sts_df[!duplicated(sts_df[, c("start", "end", "SYMBOL")]), ]
    statistiche <- subset(
        statistiche,
        select = -c(width, strand, seqnames.1, start.1, end.1, width.1, strand.1)
    )

    colnames(statistiche)[1:3] <- c("chromosome", "start", "end")
    
    # Merge with original gene list
    merge_g <- dplyr::full_join(target_bed, statistiche, by = "SYMBOL")

    # Get sample column names and ensure numeric
    col_name <- colnames(merge_g)
    col.sub <- col_name[grepl("sample_", col_name)]
    merge_g[col.sub] <- sapply(merge_g[col.sub], as.numeric)

    # CRITICAL: Calculate interval widths (number of bases in each coverage block)
    # For 1-based coordinates: width = end - start + 1
    merge_g$width <- merge_g$end.x - merge_g$start.x + 1

    # Prepare data for statistics per sample
    x <- lapply(col.sub, function(i) {
        x_df <- merge_g %>%
            dplyr::select("chromosome", "start.x", "end.x", "SYMBOL", "width", tidyselect::all_of(i))
        colnames(x_df)[6] <- "value"  # Rename coverage column
        x_df$sample <- i
        x_df
    })

    # COMPLETELY FIXED: Proper statistical calculation
    statistical_operation <- function(df) {
        # Remove rows with NA values
        df <- df %>% dplyr::filter(!is.na(value), !is.na(width))
        
        if (nrow(df) == 0) {
            return(data.frame(
                SYMBOL = character(),
                sample = character(),
                Total_bases = numeric(),
                Mean_coverage = numeric(),
                Median_coverage = numeric(),
                number_of_intervals_under_20x = numeric(),
                bases_under_20x = numeric(),
                percentage_bases_under_20x = numeric(),
                stringsAsFactors = FALSE
            ))
        }
        
        df %>%
            dplyr::group_by(SYMBOL, sample) %>%
            dplyr::summarize(
                # Total number of bases covered
                Total_bases = sum(width, na.rm = TRUE),
                
                # FIXED: Weighted mean coverage
                # Each interval contributes proportionally to its size
                # Example: [200bp @ 50x] + [100bp @ 10x] = (200*50 + 100*10)/(200+100) = 36.67x
                Mean_coverage = sum(value * width, na.rm = TRUE) / sum(width, na.rm = TRUE),
                
                # Median of interval depths (not weighted by size)
                Median_coverage = median(value, na.rm = TRUE),
                
                # Count of coverage intervals with depth < 20x
                number_of_intervals_under_20x = sum(value < 20, na.rm = TRUE),
                
                # FIXED: Total bases with coverage < 20x (weighted by interval size)
                # Example: 2 intervals under 20x with sizes 50bp and 100bp = 150 bases
                bases_under_20x = sum(width[value < 20], na.rm = TRUE),
                
                # FIXED: Percentage of BASES (not intervals!) under 20x
                # OLD WRONG FORMULA: (count_intervals / sum_of_depths) * 100
                # NEW CORRECT FORMULA: (bases_under_threshold / total_bases) * 100
                # Example: 150 bases under 20x / 500 total bases = 30%
                percentage_bases_under_20x = (sum(width[value < 20], na.rm = TRUE) / 
                                              sum(width, na.rm = TRUE)) * 100,
                
                .groups = "drop"
            ) %>%
            as.data.frame() %>%
            dplyr::mutate_if(is.numeric, round, 3)
    }

    out_r <- do.call(rbind, lapply(x, statistical_operation))

    cat("\n=== STATISTICS METHOD ===\n")
    cat("Weighted statistics accounting for interval sizes:\n\n")
    cat("Mean Coverage:\n")
    cat("  Formula: sum(depth × interval_size) / sum(interval_size)\n")
    cat("  Example: (50×200bp + 10×100bp) / 300bp = 36.67x\n\n")
    cat("Percentage Bases Under 20x:\n")
    cat("  Formula: sum(interval_size where depth<20) / total_bases × 100\n")
    cat("  Example: 100bp under threshold / 300bp total = 33.33%\n\n")
    cat("This correctly weights each coverage block by its genomic span.\n")
    
    cat(paste("\nStatistics computation took:",
              round(difftime(Sys.time(), start_time, units = "mins"), 2),
              "minutes\n\n"))

    return(out_r)
})

