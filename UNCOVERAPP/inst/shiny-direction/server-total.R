#server-preprocess.R
require(OrganismDbi)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
require(org.Hs.eg.db)
require(rlist)
require(Rsamtools)


ff=reactiveValues()
gene_list <- reactive ({
   file_gene<-input$gene1
   if (is.null(input$gene1)) {return(NULL)}
   if(input$type_file=="target_bed"){
      tmp_gene1<-read.table(input$gene1$datapath, stringsAsFactors = F)
      tmp_gene<- tmp_gene1[1:4]
      colnames(tmp_gene)<- c('chr', 'start', 'end', 'SYMBOL')
      #print(head(tmp_gene))
      ff<- tmp_gene }else{
         tmp_gene<- scan(file_gene$datapath, character(), quote = "")
         ff<- tmp_gene}
   #print(tmp_gene)
})

# list_bam<- reactive({
#    if (is.null(input$bam_list)) {return(NULL)}
#    tmp_bam <- scan(input$bam_list$datapath, character(), quote = "")
#    #list_bam= (tmp_bam)
#    #print(tmp_bam)
# })
# Lista file di coverage (bam o bed)
list_coverage <- reactive({
    if (is.null(input$bam_list)) return(NULL)
    scan(input$bam_list$datapath, character(), quote = "")
})


all_gene<- reactive({ if (input$Genome == "hg19"){
   all_gene<- TxDb.Hsapiens.UCSC.hg19.knownGene}else{
      all_gene<-TxDb.Hsapiens.UCSC.hg38.knownGene}
})

no_entrID<- reactive({
   if (is.null(gene_list()))
      return(NULL)
   if (is.null(list_coverage()))
      return(NULL)
   if(input$type_file=="target_bed"){
      errore<- data.frame(matrix(ncol = 2, nrow = 0))
   }else{# STAGE 1: Convert gene names to ENTREZ IDs using org.Hs.eg.db
    cat("Stage 1: Converting gene names to ENTREZ IDs using org.Hs.eg.db...\n")

    # Handle 1:many mappings by taking the first ENTREZ ID for each gene
    my_gene_name <- OrganismDbi::select(org.Hs.eg.db, key = gene_list(), columns = "ENTREZID", keytype = "ALIAS") %>%
      tibble::as_tibble() %>%
      dplyr::filter(!is.na(ENTREZID)) %>%
      dplyr::group_by(ALIAS) %>%
      dplyr::slice(1) %>%  # Take first ENTREZ ID if multiple exist
      dplyr::ungroup()

    # Check for genes not found in org.Hs.eg.db
    found_genes_stage1 <- unique(my_gene_name$ALIAS)
    not_found_stage1 <- base::setdiff(gene_list(), found_genes_stage1)

    if (length(not_found_stage1) > 0) {
      cat("WARNING: These genes were not found in org.Hs.eg.db and will be excluded:\n")
      for(gene in not_found_stage1) {
      cat(paste("  - Gene '", gene, "' not found in org.Hs.eg.db\n", sep=""))
    }
  cat("\n") #
}

if (nrow(my_gene_name) == 0) {
  stop("ERROR: No valid genes found in org.Hs.eg.db. Please check your gene names.")
}

cat(paste("Stage 1 result:", nrow(my_gene_name), "gene-ID mappings found from", length(found_genes_stage1), "unique genes\n\n"))

# STAGE 2: Check if ENTREZ IDs exist in genome reference database
cat("Stage 2: Validating genes against reference genome database...\n")

all_entrez_ids <- unique(my_gene_name$ENTREZID)
genome_keys <- AnnotationDbi::keys(all_gene(), "GENEID")  # Fixed: added () to call the reactive

# Find genes present in both databases
valid_entrez_ids <- base::intersect(all_entrez_ids, genome_keys)
invalid_entrez_ids <- base::setdiff(all_entrez_ids, genome_keys)

if (length(invalid_entrez_ids) > 0) {
  # Get gene names for invalid ENTREZ IDs
  invalid_genes_stage2 <- my_gene_name %>% 
    dplyr::filter(ENTREZID %in% invalid_entrez_ids) %>%
    dplyr::pull(ALIAS) %>%
    unique()
  
  cat("WARNING: These genes were found in org.Hs.eg.db but not in the reference genome database and will be excluded:\n")
  for(gene in invalid_genes_stage2) {
    cat(paste("  - Gene '", gene, "' not found in reference genome database (", input$Genome, ")\n", sep=""))  # Fixed: changed 'genome' to 'input$Genome'
  }
  cat("\n")
} else {
  invalid_genes_stage2 <- character(0)
}

# Filter to keep only valid genes
my_gene_name <- my_gene_name %>%
  dplyr::filter(ENTREZID %in% valid_entrez_ids)

ID <- valid_entrez_ids

if (length(ID) == 0) {
  stop("ERROR: No genes passed both validation stages. Please check your gene names and genome version.")
}

cat(paste("Stage 2 result:", length(valid_entrez_ids), "genes validated against reference genome\n"))
cat(paste("FINAL RESULT: Proceeding with", length(ID), "validated genes for analysis\n"))

# Print summary of excluded genes
total_excluded <- length(not_found_stage1) + length(invalid_genes_stage2)
if (total_excluded > 0) {
  cat(paste("SUMMARY: ", total_excluded, " genes excluded (", 
            length(not_found_stage1), " from stage 1, ", 
            length(invalid_genes_stage2), " from stage 2)\n", sep=""))
}
cat("=== GENE VALIDATION COMPLETE ===\n\n")

# Return error data frame if there are issues, otherwise return empty data frame
if (length(ID) == 0) {
  errore <- data.frame(
    Gene = c("No valid genes found"),
    Issue = c("All genes failed validation"),
    stringsAsFactors = FALSE
  )
  return(errore)
} else {
  errore <- data.frame(matrix(ncol = 2, nrow = 0))
  return(errore)
}
  }
})



#       gene_list1<- gene_list()
#       my_gene_name<- OrganismDbi::select(org.Hs.eg.db, key= gene_list1, columns="ENTREZID",
#                                          keytype="ALIAS")
#       #print(my_gene_name)
#       ID<- my_gene_name$ENTREZID
#       b<- c()
#       for (i in ID){
#          if( ! is.element(i, keys(all_gene(), "GENEID")))
#             b[i]<- i
#       }
#       errore<- subset(my_gene_name, my_gene_name$ENTREZID %in% b)
#       colnames(errore)[1]= " the following HGNC official gene names are
#     unrecognizable, please choose a correct name and reload a file"
#       return(errore)
#       #print(head(errore,n=12))
#    }
# })


for_bed<- reactive({
   if (is.null(gene_list()))
      return(NULL)
   if (is.null(list_coverage()))
      return(NULL)
   if(nrow(no_entrID())!=0)
      return(no_entrID())
   if(input$type_file=="target_bed"){
      for_bed<- gene_list() }else{
         gene_list1<- gene_list()
         
         # Get valid gene mappings (same logic as in no_entrID but simplified)
         my_gene_name <- OrganismDbi::select(org.Hs.eg.db, key = gene_list1, columns = "ENTREZID", keytype = "ALIAS") %>%
           tibble::as_tibble() %>%
           dplyr::filter(!is.na(ENTREZID)) %>%
           dplyr::group_by(ALIAS) %>%
           dplyr::slice(1) %>%  # Take first ENTREZ ID if multiple exist
           dplyr::ungroup()
         
         # Get valid ENTREZ IDs that exist in the genome database
         all_entrez_ids <- unique(my_gene_name$ENTREZID)
         genome_keys <- AnnotationDbi::keys(all_gene(), "GENEID")
         valid_entrez_ids <- base::intersect(all_entrez_ids, genome_keys)
         
         # Filter to keep only valid genes
         my_gene_name <- my_gene_name %>%
           dplyr::filter(ENTREZID %in% valid_entrez_ids)
         
         ID <- valid_entrez_ids
         pre<- data.frame()
         
         for (i in ID){
            # Check if the ENTREZ ID is valid before using it
            if (i %in% genome_keys) {
               txid <- OrganismDbi::select(all_gene(), i, "TXNAME", "GENEID")[["TXNAME"]]
               cds <-OrganismDbi::exonsBy(all_gene(), by="tx", use.names=TRUE)
               coor<- as.data.frame(cds[names(cds) %in% txid])
               coordinate<- cbind(coor,i)
               colnames(coordinate)[11]<- "ENTREZID"
               cood<- dplyr::inner_join(coordinate, my_gene_name, by="ENTREZID")
               pre<- rbind(cood,pre)
               pre$start<- pre$start -10
               pre$end<- pre$end +10
            }
         }
         #print(head(pre))
         for_bed<- data.frame()
         for (row in 1:nrow(pre)){
            sel_cc<- data.frame(as.character(pre$seqnames[row]),
                                pre$start[row], pre$end[row],pre$ALIAS[row],
                                stringsAsFactors = FALSE)
            for_bed<- rbind(for_bed, sel_cc)

         }
         colnames(for_bed)<- c('chr', 'start', 'end', 'SYMBOL')
         for_bed<- unique(for_bed)
         if(input$notation == "number"){
            for_bed$chr<- gsub("^.{0,3}", "", for_bed$chr, perl =TRUE)}
            
         return(for_bed)
         print(for_bed)
      }
})


coverage_input <- reactive({
    if (is.null(gene_list())) return(NULL)
    if (!is.null(no_entrID()) && nrow(no_entrID()) != 0) return(no_entrID())
    if (is.null(list_coverage())) return(NULL)

    if (input$type_coverage == "bam") {
        # --- calcolo coverage dai BAM ---
        for_grange <- GenomicRanges::makeGRangesFromDataFrame(for_bed(), keep.extra.columns = TRUE)
        print("for_grange") 
        print(for_grange)
        param <- Rsamtools::ScanBamParam(which=for_grange)
        paste("param")
        print(param)

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
        # Normalize chromosome notation once (canonical: chr*)
        pp$seqnames <- paste0("chr", sub("^chr", "", as.character(pp$seqnames)))
        g<-pp
         print("g")
        print(g)
        return(pp)
    }
    if (input$type_coverage == "bed") {
      # 1. Get target gene regions (same as BAM)
       for_grange <- GenomicRanges::makeGRangesFromDataFrame(
          for_bed(), 
           keep.extra.columns = TRUE
      )
    
    # 2. Read and filter each BED file
       df_list <- list()
      for (file in list_coverage()) {
        
          cat("\n=== READING BED FILE ===\n")
         cat(paste("File:", file, "\n"))
        
        # Read BED - assuming 4 columns: chr, start, end, coverage
         # Don't assume header exists in BED files
         df <- read.table(file, 
                         header = FALSE,
                         stringsAsFactors = FALSE,
                         colClasses = c("character", "integer", "integer", "numeric"))
        
        # Assign column names
           colnames(df) <- c("seqnames", "start", "end", "count")

          cat("Before filtering - rows:", nrow(df), "\n")
           print(head(df))
        
        # Convert to GRanges for filtering
         df_gr <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
         # Right after creating df_gr and before filtering
         cat("\n=== COORDINATE COMPARISON ===\n")
         cat("BED Coverage Ranges:\n")
         cat("  Chromosome:", as.character(unique(seqnames(df_gr))), "\n")
         cat("  Min position:", min(start(df_gr)), "\n")
         cat("  Max position:", max(end(df_gr)), "\n")
         cat("  Example regions:\n")
         print(head(df_gr, 3))

         cat("\nGene Target Ranges:\n")
         cat("  Chromosome:", as.character(unique(seqnames(for_grange))), "\n")
         cat("  Min position:", min(start(for_grange)), "\n")
         cat("  Max position:", max(end(for_grange)), "\n")
         cat("  Example regions:\n")
         print(head(for_grange, 3))

         cat("\n=== CHECKING CHROMOSOME NAMES ===\n")
         cat("BED seqlevels:", seqlevels(df_gr), "\n")
         cat("Gene seqlevels:", seqlevels(for_grange), "\n")
         cat("Do they match?", any(seqlevels(df_gr) %in% seqlevels(for_grange)), "\n")
         cat("============================\n\n")

         # Now do the filtering
         df_gr_filtered <- subsetByOverlaps(df_gr, for_grange, type="any")

        
           # Filter to only gene regions (like ScanBamParam does for BAM)
         df_gr_filtered <- IRanges::subsetByOverlaps(df_gr, for_grange, type="any")
        
          cat("After filtering - regions:", length(df_gr_filtered), "\n")
        
          # Convert back to data frame
           df_filtered <- as.data.frame(df_gr_filtered)
        
        # Keep only the columns we need: seqnames, start, end, count
        df_filtered <- df_filtered[, c("seqnames", "start", "end", "count")]
        
         cat("Filtered data frame structure:\n")
         print(str(df_filtered))
         print(head(df_filtered))
         cat("===========================\n\n")
        
           df_list <- append(df_list, list(df_filtered))
       }
    
    # 3. Merge all samples (exactly like BAM processing!)
       pp <- Reduce(function(...) merge(..., by=c("seqnames","start","end"), all=TRUE), df_list)
     pp[is.na(pp)] <- 0
    
    # 4. CRITICAL: Convert types to match BAM output exactly
    # BAM produces: seqnames=factor, start=int, end=int, count=int
       pp$seqnames <- factor(pp$seqnames)
       pp$start <- as.integer(pp$start)
       pp$end <- as.integer(pp$end)
    
    # Convert all count columns to integer
       count_cols <- grep("^count", colnames(pp), value = TRUE)
     pp[count_cols] <- lapply(pp[count_cols], as.integer)
    
    # Verify output structure matches BAM
       cat("\n=== FINAL OUTPUT STRUCTURE ===\n")
       cat("Column types:\n")
       print(sapply(pp, class))
       cat("First rows:\n")
       print(head(pp))
       cat("===============================\n\n")
    
      # Normalize chromosome notation once (canonical: chr*)
      pp$seqnames <- paste0("chr", sub("^chr", "", as.character(pp$seqnames)))
      g <- pp
      print("g")
      print(g)
    
      return(pp)
   }



        
    
})

# pileup_input<- reactive({
#    if (is.null(gene_list()))
#       return(NULL)
#    if (is.null(list_bam()))
#       return(NULL)
#    if(nrow(no_entrID())!=0)
#       return(no_entrID())

#    for_grange<-GenomicRanges::makeGRangesFromDataFrame(for_bed(),
#                                                        keep.extra.columns = TRUE)
#    #print(head(for_grange))
#    param <- Rsamtools::ScanBamParam(which= for_grange)

#    p_param <- Rsamtools::PileupParam(distinguish_nucleotides=FALSE,
#                                      distinguish_strands=FALSE,
#                                      min_mapq=as.numeric(input$MAPQ),
#                                      min_base_quality=as.numeric(input$base_qual),
#                                      min_nucleotide_depth=1,
#                                      max_depth=150000)
#    df= list()
#    for (i in list_bam()){
#       pileup_df= Rsamtools::pileup(i, scanBamParam=param, pileupParam=p_param)

#       df=rlist::list.append(df, pileup_df)
#    }
#    lst1 <- lapply(df, function(x) transform(x[,-5]))
#    lst2<- lapply(lst1, function(x) transform(x[!duplicated(x),]))

#    riarrange.df <- function(list_df){
#       require(dplyr)
#       list_df %>%
#          dplyr::mutate(end= pos) %>%
#          dplyr::group_by(seqnames,pos,end) %>%
#          dplyr::summarise(count=sum(count)) %>%
#          dplyr::arrange() }

#          #dplyr::summarise(value= as.numeric(paste(sum(count))),
#                           #counts=paste(count, sep=':',collapse=';'
   
#    lst3<- lapply(lst2, riarrange.df)

#    pp<- Reduce(function(...) merge(...,by= c("seqnames", "pos", "end")), lst3)
#    pp[is.na(pp)] <- 0
#    pp<- as.data.frame(pp)
#    if(input$notation == "number"){
#       for (i in pp[1]){
#          Chromosome<-paste("chr", i, sep="")
#          pp<- cbind(Chromosome, pp)
#          pp[,2]<-NULL
#       }
#    }
#    colnames(pp)<-NULL
#    return(pp)
# })


name_sample<- reactive({
   if (is.null(list_coverage())) return(character(0))
   tools::file_path_sans_ext(basename(list_coverage()))
})

stat_summ<- reactive({
   if (is.null(gene_list()))
      return(NULL)
   if (is.null(list_coverage()))
      return(NULL)
   if(nrow(no_entrID())!=0)
      return(no_entrID())
   ppinp<- as.data.frame(coverage_input())
   print(ppinp)
   colnames(ppinp)[1:3]<-c("seqnames","start","end")
  n<- length(colnames(ppinp)[-1:-3])
  samples <- name_sample()
  colnames(ppinp)[-1:-3] <- paste0("sample_", samples)[seq_len(n)]

  #colnames(ppinp)[-1:-2]<-paste0(c("sample_"), samples[1:n])
   for_grange<-GenomicRanges::makeGRangesFromDataFrame(for_bed(),
                                                       keep.extra.columns = TRUE)

   for_range_pp=GenomicRanges::makeGRangesFromDataFrame(ppinp,
                                                        keep.extra.columns = TRUE)
   tp<- GenomicRanges::findOverlaps(query= for_range_pp,
                                    subject = for_grange, type="any",
                                    select = "all")
   sts_df <- data.frame(for_range_pp[queryHits(tp),], for_grange[subjectHits(tp),])
   print(sts_df)
   ### Code diagnostics 7/10/2025
   # Dopo: sts_df <- data.frame(for_range_pp[queryHits(tp),], for_grange[subjectHits(tp),])
   cat("\n=== OVERLAP DIAGNOSTICS ===\n")
   cat(paste("Number of overlaps found:", length(tp), "\n"))
   cat(paste("Unique genes with coverage:", length(unique(sts_df$SYMBOL)), "\n"))
   cat("Genes WITH coverage:\n")
   print(unique(sts_df$SYMBOL))
   cat("Genes WITHOUT coverage:\n")
   genes_no_coverage <- setdiff(unique(for_bed()$SYMBOL), unique(sts_df$SYMBOL))
   print(genes_no_coverage)
   cat("Coverage chromosomes:\n")
   print(unique(ppinp$seqnames))
   cat("Gene chromosomes:\n")
   print(table(for_bed()$chr))
   cat("===========================\n\n")
   ### end

   statistiche<- sts_df[!duplicated(sts_df$start),]
   statistiche<- subset(statistiche,
                        select = -c(width, strand, seqnames.1,
                                    start.1, end.1, width.1, strand.1))

   colnames(statistiche)[1:3]<- c("chromosome","start","end")
   merge_g<- dplyr::full_join(for_bed(),statistiche, by="SYMBOL")
   #merge_g <- dplyr::full_join(for_bed(), statistiche, by="SYMBOL")
   
   #col.sub<- colnames(merge_g[, grepl("sample_" , names(merge_g))])
   col_name= colnames(merge_g)

   col.sub= col_name[grepl("sample_", col_name)]
  merge_g[col.sub] <- sapply(merge_g[col.sub],as.numeric)

   x<- list()
   for (i in col.sub){
      x_df<- merge_g %>%
      dplyr::select("chromosome","start.x","end.x", "SYMBOL",i)
      #print(x_df)
      colnames(x_df)[5]<- "value"
      x_df$sample<- paste0(i)
      x<- list.append(x, x_df)
   }
   statistical_operation= function(df){
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

   out_r<- do.call(rbind, lapply(x, function(x) statistical_operation(x)))
   return(out_r)


})

#server-annotation.R
filtered_low_nucl<- reactive ({
  if (is.null(mysample()))
    return(NULL)
  
  df <- mysample()
  thr <- input$coverage_co
  
  # Get chromosome value - use input directly to avoid reactive timing issues
  chr_val <- input$Chromosome
  
  cat("Filtering for chromosome:", chr_val, "with coverage threshold:", thr, "\n")
  cat("Sample data dimensions:", dim(df), "\n")
  cat("Available chromosomes in data:", unique(df$chromosome), "\n")
  
  if (identical(thr, "all")) {
    result <- dplyr::filter(df, chromosome == chr_val)
  } else {
    result <- dplyr::filter(df, chromosome == chr_val,
                  coverage <= as.numeric(thr))
  }
  
  cat("Filtered result dimensions:", dim(result), "\n")
  if (nrow(result) > 0) {
    cat("Sample of filtered data:\n")
    print(head(result))
  } else {
    cat("No data found for chromosome", chr_val, "with coverage <=", thr, "\n")
  }
  
  return(result)
})

intBED<- reactive({
  cat("=== STARTING intBED FUNCTION ===\n")
  
  if (is.null(filtered_low_nucl())) {
    cat("ERROR: filtered_low_nucl() returned NULL\n")
    return(NULL)
  }
  
  bedA<- filtered_low_nucl()
  cat("bedA (coverage data) dimensions:", dim(bedA), "\n")
  cat("bedA sample data:\n")
  print(head(bedA))

  # TEST OVERRIDE: prefer custom local files or env vars for annotation during testing
  test_hg19 <- "/home/anna/uncoverappLib/traials/geneApp/R/sorted_hg19.bed.gz"
  test_hg38 <- "/home/anna/uncoverappLib/traials/geneApp/R/sorted_hg38.bed.gz"
  env_hg19 <- Sys.getenv("UNCOVERAPP_HG19_ANNOTATION", unset = "")
  env_hg38 <- Sys.getenv("UNCOVERAPP_HG38_ANNOTATION", unset = "")

  file.name <- NULL
  # 1) Environment variable overrides (if set)
  if (identical(input$UCSC_Genome, "hg19") && nzchar(env_hg19) && file.exists(env_hg19)) {
    file.name <- env_hg19
  }
  if (identical(input$UCSC_Genome, "hg38") && nzchar(env_hg38) && file.exists(env_hg38)) {
    file.name <- env_hg38
  }
  # 2) Test files in repo (if present)
  if (is.null(file.name) && identical(input$UCSC_Genome, "hg19") && file.exists(test_hg19)) {
    file.name <- test_hg19
  }
  if (is.null(file.name) && identical(input$UCSC_Genome, "hg38") && file.exists(test_hg38)) {
    file.name <- test_hg38
  }
  # 3) Fallback to default cache selection
  if (is.null(file.name)) {
    m <- uncoverappLib::getAnnotationFiles()
    # Separate data files (.gz) from index files (.tbi)
    data_files <- m[grepl("\\.gz$", m) & !grepl("\\.tbi$", m)]
    # Selection rules:
    # - If multiple data files exist (e.g., hg19 + hg38), pick by genome string.
    # - Else, use the only data file present (ignore genome choice).
    if (length(data_files) > 1) {
      idx <- if (input$UCSC_Genome == "hg19") {
        grep("hg19|GRCh37", data_files, ignore.case = TRUE)
      } else {
        grep("hg38|GRCh38", data_files, ignore.case = TRUE)
      }
      if (length(idx) == 0) idx <- 1
      file.name <- data_files[idx[1]]
    } else if (length(data_files) == 1) {
      file.name <- data_files[1]
    } else {
      # Fallback: if nothing matched (unexpected), try the first returned path
      file.name <- m[1]
    }
  }
  # Print to console which file is used (useful during testing)
  cat("Using annotation file:", file.name, "\n")

  #second and tirth columns are hg19 positions
  cat("\n=== PROCESSING QUERY ===\n")
  query <- c(input$query_Database)
  cat("Original query string:", query, "\n")
  
  query.regions= read.table(text=gsub("[:-]+", " ", query, perl=TRUE),
             header=FALSE, col.names = c("chr", "start", "end"))
  
  if (is.null(query.regions)) {
    cat("ERROR: Failed to parse query regions\n")
    return(NULL)
  }
  
  cat("Parsed query regions:\n")
  print(query.regions)
  cat("\n=== RUNNING TABIX QUERY ===\n")
  result<- try({
    cat("Creating GRanges from query regions...\n")
    fq= GenomicRanges::makeGRangesFromDataFrame(query.regions, keep.extra.columns = TRUE)
    cat("GRanges created successfully\n")
    
    cat("Running scanTabix with file:", file.name, "\n")
    res <- Rsamtools::scanTabix(file.name, param=fq)
    cat("Tabix completed\n")
    
    lengths_result <- sapply(res, length)
    cat("Tabix result lengths:", lengths_result, "\n")
    dff <- Map(function(elt) {
      if (length(elt) == 0) {
        return(data.frame())
      }
      read.csv(textConnection(elt), sep="\t", header=FALSE, stringsAsFactors = FALSE)
    }, res)
    
    # Check if we have any data
    valid_dfs <- dff[sapply(dff, nrow) > 0]
    if (length(valid_dfs) == 0) {
      cat("No annotation data found for query region\n")
      return(NULL)
    }
    
    bedB <- do.call(rbind, valid_dfs)

  })
  if ("try-error" %in% class(result)) {
    err_msg <- 'no coordinates recognized'
    cat("Error in annotation processing:", err_msg, "\n")
    return(NULL)
  }
  
  # Check if bedB exists and has data
  if (is.null(bedB) || nrow(bedB) == 0) {
    cat("No annotation data to process\n")
    return(NULL)
  }

  print(head(bedB))
  # Check actual number of columns and set names accordingly
  ncols <- ncol(bedB)
  cat("Number of columns in annotation data:", ncols, "\n")
  
  if (ncols == 19) {
    colnames(bedB)<- c ('Chromo', 'start','end','REF','ALT',
                        'dbsnp','GENENAME', 'PROTEIN_ensembl', 'field9',
                        'MutationAssessor','SIFT','Polyphen2',
                        'M_CAP','CADD_PHED','AF_gnomAD','ClinVar',
                        'clinvar_MedGen_id','HGVSc_VEP','HGVSp_VEP')
  } else {
    # Fallback: use generic column names
    colnames(bedB) <- paste0("V", 1:ncols)
  }
  str(bedB)
  
  # For this data structure, start and end are already in the right columns
  # No need to rename - they're already named 'start' and 'end'
  cat("Column names after assignment:", colnames(bedB), "\n")

  for (i in bedB[1]){
    Chromosome<-paste("chr", i, sep="")
    bedB<- cbind(Chromosome, bedB)
    bedB[,2]<-NULL
  }
  bedB$Chromosome= as.character(bedB$Chromosome)
  bedB$AF_gnomAD= suppressWarnings(as.numeric(bedB$AF_gnomAD))
  bedB$CADD_PHED= suppressWarnings(as.numeric(bedB$CADD_PHED))
  intersectBedFiles.GR <- function(bed1,bed2) {
    require(GenomicRanges)
    require(IRanges)
    bed1 <- makeGRangesFromDataFrame(bedA,ignore.strand = TRUE,
                                     keep.extra.columns = TRUE)
    bed2 <- makeGRangesFromDataFrame(bedB,ignore.strand = TRUE,
                                     keep.extra.columns = TRUE)
    tp<- findOverlaps(query = bed1, subject = bed2, type="any")
    intersect_df = data.frame(bed1[queryHits(tp),], bed2[subjectHits(tp),])
    return(intersect_df)
  }
  intersect_df<- intersectBedFiles.GR(bedA, bedB)
  return(intersect_df)
})

#make reactive dataframe

condform_table<- reactive ({
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "table construction in progress",
               detail = 'This may take a while', value = 0)
  Sys.sleep(0.1)
  library(condformat)
  if (is.null(intBED()))
    return(NULL)
  condformat(intBED()) %>%
    rule_fill_discrete(ClinVar, expression= ClinVar !=".",
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(CADD_PHED, expression= CADD_PHED > 20,
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(MutationAssessor, expression=  MutationAssessor =='H',
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(M_CAP, expression=  M_CAP =='D',
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(AF_gnomAD, expression=
                         ifelse(is.na(AF_gnomAD) | AF_gnomAD < 0.5,
                                'TRUE','FALSE') ,
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(c(start, end),
                       expression = grepl("H|M", MutationAssessor) &
                         ClinVar !="." &  AF_gnomAD < 0.5 ,
                       colours= c("TRUE"= "yellow", "FALSE"= ""))%>%
    rule_css(c(start, end),
             expression = ifelse(grepl("H|M", MutationAssessor) &
                                   ClinVar !="." &  AF_gnomAD < 0.5,
                                 "red", "green"),
             css_field = "color")
})
#server-plots.R

#strack<- reactive({if (input$UCSC_Genome == "hg19"){
#  Gviz::SequenceTrack(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, chromosome = Chromosome())}
#  else{
#    Gviz::SequenceTrack(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, chromosome = Chromosome())}
#})

dtrack1<- reactive ({grcoverage <- filtered_low()
dt1<-Gviz::DataTrack(range= grcoverage,type = "histogram", fill.histogram =  "red", col.histogram="NA",
               genome = input$UCSC_Genome,name = "Seq. Depth")
})

dtrackHigh<- reactive ({ grcoverage_high <- filtered_high()
dt2<- Gviz::DataTrack(range = grcoverage_high,type = "histogram", fill.histogram = "dodgerblue",
                col.histogram="NA", genome = input$UCSC_Genome,name = "Seq. Depth")
})

itrack <- reactive ({
  i= Gviz::IdeogramTrack(genome =input$UCSC_Genome , chromosome = Chromosome())
})

grtrack  <-reactive({
  ggT<-Gviz::GeneRegionTrack (txdb(), chromosome  = Chromosome() ,  start  =  St() ,  end  =  En(),
                        showId  =  TRUE , transcriptAnnotation="symbol",
                        name  =  " Gene Annotation ")
  return(ggT)
})
#Prepration for all gene coverage plot
p1<- reactive ({
  start_gene= coord()$start[coord()$seqnames ==Chromosome()]
  if (length(start_gene) > 1)
    start_gene= start_gene[1]
  print(start_gene)
  end_gene= coord()$end[coord()$seqnames ==Chromosome()]
  if (length(end_gene) > 1)
    end_gene= end_gene[1]
  print(end_gene)

  ot <- Gviz::OverlayTrack(trackList = list(dtrack1(), dtrackHigh()))
  gtrack <- Gviz::GenomeAxisTrack()
  #ylims <- extendrange(range(c(values(dtrack1()), values(dtrackHigh()))))
  ylims <- grDevices::extendrange(range(c(dtrack1()@data), dtrackHigh()@data))
  gr_ex_track <- Gviz::GeneRegionTrack(txdb(),
                                       chromosome = Chromosome(),
                                       start = start_gene, end = end_gene,
                                 showId = TRUE,
                                 name = "Gene Annotation")
  p1= Gviz::plotTracks(list(itrack(), gtrack, ot, gr_ex_track),
                 from = start_gene,
                 to=end_gene,reverseStrand = TRUE,
                 ylim=ylims,type = "histogram",baseline=input$coverage_co,
                 col.baseline = "black",lwd.baseline=0.3,
                 extend.left = 0.5, extend.right = 200)

})

#table of uncovered exons

table1<- reactive({
  if (is.null(filtered_low()))
    return(NULL)
  f.low=filtered_low()
  f.low[,'new']='NA'
  f.low$new<- ifelse(sapply(f.low$start, function(p)
    any(exon_gp()$start <= p & exon_gp()$end >= p)), "YES", "out")
  l.coverage= as.data.frame(f.low[f.low$new=='YES',])
  validate(
    need(nrow(l.coverage) >0, "ALL EXONS ARE COVERED
         UNDER YOUR CHOOSE THRESHOLD"))

  x= l.coverage$start
  getValue3 <- function(x, data) {
    tmp <- data %>%
      dplyr::filter(start <= x, x <= end) %>%
      dplyr::filter(number_of_transcript == input$transcript_id)
    return(tmp$exon_rank)
  }

  a_exon=sapply(x, getValue3, data=exon_gp())
  exon=unlist(lapply(a_exon, function (x) ifelse(length (x) > 0, x, NA)))
  exon.df= as.data.frame(exon)
  if (is.null(exon.df))
    return(NULL)
  df.l= cbind(l.coverage, exon.df)
  t1=table(df.l$exon)
  df.t1= as.data.frame(t1)
  colnames(df.t1)= c('exon','uncovered positions')
  return(df.t1)
})

#output$df.l<- renderDataTable({
 # table1()
#})

#Preparation of exon zooom plot

p2<- eventReactive(input$button5, {
  gname =input$Gene_name
  if (is.null(gname))
    return(NULL)
  disable("button5")
  shinyjs::show("text1")
  Sys.sleep(0.1)
  gtrack <- Gviz::GenomeAxisTrack()
  ot <- Gviz::OverlayTrack(trackList = list(dtrack1(), dtrackHigh()))
  ylims <- grDevices::extendrange(range(c(dtrack1()@data), dtrackHigh()@data))
  one_trascript= grtrack()[grtrack()@range$transcript== id()]
  print(head(one_trascript))
  grtrack_symbol <- Gviz::GeneRegionTrack(one_trascript@range,
                                          chromosome = Chromosome(),
                                    start = St(),
                                    end = En(),
                                    showId = TRUE, exonAnnotation="exon",
                                    name = "Gene Annotation & Symbol")
  grtrack_range <- grtrack_symbol@range
  range_mapping <- OrganismDbi::select(Homo.sapiens,
                                       keys = mcols(grtrack_range)$symbol,
                                       keytype = "TXNAME",
                                       columns = c("ENTREZID", "SYMBOL"))
  library(stringr)
  new_symbols <- with(range_mapping,str_c(SYMBOL, " (", TXNAME, ")", sep = ""))
  symbol(grtrack_symbol) <- new_symbols
  shinyjs::enable("button5")
  shinyjs::hide("text1")
  p2=Gviz::plotTracks(list(itrack(), gtrack, ot, grtrack_symbol),
                      ylim=ylims,xlim=NULL,exonAnnotation="exon1",
                from = St(), to = En(), reverseStrand = TRUE,
                background.panel = "#FFFEDB", background.title = "darkblue",
                baseline=input$coverage_co,
                col.baseline = "black",lwd.baseline=0.3,
                extend.left = 0.5, extend.right = 100,
                main= paste0('exon',input$exon_number))

})

observeEvent(input$button5,{
  output$ens<- renderPlot({
    validate(
      need(ncol(mydata()) != "0", "Unrecognized data set:
           Please load your file"))
    options(shiny.sanitize.errors = TRUE)
    progress <- shiny::Progress$new()
    for (i in 1:40) {
      progress$set(message = "Please wait a few minutes:
                   Making plot", detail = 'This may take a while', value = i)
      Sys.sleep(1.0)
    }
    on.exit(progress$close())
    Sys.sleep(0.1)
    p2()
  })
})

observeEvent(input$remove5,{
  shinyjs::hide(p2())
  output$ens<- NULL
})

observeEvent(input$btn_run,{
  Sys.sleep(5)
})

####output plot with reference sequqnce, few bases##########

p3<- reactive({
  ot <- Gviz::OverlayTrack(trackList = list(dtrack1(), dtrackHigh()))
  gtrack <- Gviz::GenomeAxisTrack()
  ylims <- grDevices::extendrange(range(c(dtrack1()@data), dtrackHigh()@data))
  #strack<- SequenceTrack(homo, chromosome = Chromosome())
  #p3=Gviz::plotTracks(list(itrack(), gtrack, ot,grtrack(), strack()),
  p3=Gviz::plotTracks(list(itrack(), gtrack, ot,grtrack()),
                from = as.numeric(input$Start_genomicPosition),
                to=as.numeric(input$end_genomicPosition),
                reverseStrand = TRUE, cex = 0.8, ylim= ylims,
                baseline=input$coverage_co,
                col.baseline = "black",lwd.baseline=0.5,
                extend.right = 100, extend.left = 0.5)
})



#server-tables.R
#Tables
tryObserve <- function(x) {
  x <- substitute(x)
  env <- parent.frame()
  observe({
    tryCatch(eval(x, env),
             error = function(e) {
               #showNotification(paste("Error: ", e$message), type = "error")
             })
  })
}

coord= eventReactive(input$button1,{
  disable("button1")
  shinyjs::show("text1")
  Sys.sleep(0.1)
  my_gene_name=OrganismDbi::select(org.Hs.eg.db,
                                   key= input$Gene_name,
                                   columns=c("ENTREZID","GENENAME", "ENSEMBL"),
                                   keytype="ALIAS")
  ID=my_gene_name$ENTREZID
  if (is.null(ID))
    return(NULL)
  all_gene= data.frame(genes(txdb()))
  pre= do.call(rbind, lapply(ID, function(x) data.frame(
    subset(all_gene, all_gene$gene_id == x), stringsAsFactors = FALSE)))
  colnames(pre)[6]= 'ENTREZID'
  info= merge(pre, my_gene_name)
  shinyjs::enable("button1")
  shinyjs::hide("text1")
  return(info)
})


observeEvent(input$button1, {
  output$ccg <-DT::renderDataTable({
    progress <- shiny::Progress$new()
    progress$set(message = "Running", detail = 'This may take a while')
    on.exit(progress$close())
    validate(
      need(input$Gene_name != "",
        "Unrecognized gene name: Please select HGNC gene name \n Click apply"))
    HGNC_org <- keys(org.Hs.eg.db, keytype = "SYMBOL")
    validate(
      need(HGNC_org[HGNC_org %in% input$Gene_name],
           "incorrect HGNC Gene symbol"))
    coord()

  })
})

tryObserve({
  if (is.null(input$button1) )
    return()
  if (is.null(coord()))
    return(NULL)
  x= as.data.frame(coord())
  Chrom= as.character(x$seqnames[1])
  updateTextInput(session, "Chromosome", value = Chrom)
  })

Chromosome<- reactive({
  xc= as.character(input$Chromosome)
  return(xc)
})



observeEvent(input$remove,{
  shinyjs::hide(coord())
  output$ccg<- NULL
})


exon_gp<-eventReactive(input$button1,{
  #require(rtracklayer)
  #ucsc <- browserSession()
  gname =input$Gene_name
  if (is.null(gname))
    return(NULL)
  else if (gname== "")
    return(NULL)
  if (input$UCSC_Genome == "hg19"){
    edb = EnsDb.Hsapiens.v75}
  else{
    edb= EnsDb.Hsapiens.v86}
  eid <- OrganismDbi::select(org.Hs.eg.db,gname,
                             "ENTREZID", "SYMBOL")[["ENTREZID"]]
  txid <- OrganismDbi::select(txdb(), eid, "TXNAME", "GENEID")[["TXNAME"]]
  cds <- cdsBy(txdb(), by="tx", use.names=TRUE)
  exoncds <- cds[names(cds) %in% txid]
  exon_Id<- as.data.frame(exoncds)
  exon_table=exon_Id  %>%
    dplyr::select(1:6,8,10)
  #exon_table=exon_Id[c(1:6,8,10)]
  colnames(exon_table)=c("number_of_transcript",
                         "type_of_transcript", "chrom", "start","end",
                         "length_of_exon", "cds_id", "exon_rank")
  return(exon_table)
  print(head(exon_table))
})

observeEvent(input$button1, {
  output$exon_pos<- DT::renderDataTable({
    options(shiny.sanitize.errors = TRUE)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "table construction in progress",
                 detail = 'This may take a while', value = 0)
    Sys.sleep(0.1)

    HGNC_org <- keys(org.Hs.eg.db, keytype = "SYMBOL")
    validate(
      need(HGNC_org[HGNC_org %in% input$Gene_name],
           "incorrect HGNC Gene symbol"))
    exon_gp()
  })
})



observeEvent(input$remove,{
  output$exon_pos<- DT::renderDataTable({
    shinyjs::js$reset(exon_gp()) })
})


tryObserve({
  if (is.null(input$exon_number) )
    return()
  message("load inputs")
  #observe({
  # if (is.null(input$exon_number))
  #  return(NULL)
  if (is.null(exon_gp()))
    return(NULL)
  exon_df= as.data.frame(exon_gp())
  one=subset(exon_df,
             exon_df$exon_rank == input$exon_number &
               exon_df$number_of_transcript == input$transcript_id)
  print(one)
  start_exon= as.numeric(one$start)
  end_exon= as.numeric(one$end)
  chr= as.numeric(gsub("\\D", "", one$chrom, perl =TRUE))
  trascript_name= as.character(one$type_of_transcript)
  updateTextInput(session, "id_t", value = trascript_name)
  updateTextInput(session, "Start_genomicPosition", value = start_exon)
  updateTextInput(session, "end_genomicPosition", value = end_exon)
  updateTextInput(session, "query_Database",
                  value = paste0(chr,":",start_exon,"-",end_exon))
})

#I define reactive START and END position

St<- reactive({
  sta= as.numeric(input$Start_genomicPosition)
  return(sta)
})

En<- reactive({
  en= as.numeric(input$end_genomicPosition)
  return(en)
})

id<- reactive({
  id_ucsc= as.character(input$id_t)
  return(id_ucsc)
})

coordinate<- reactive({
  coo=paste0(chr,":",start_exon,"-",end_exon)
  return(coo)
})
#server-maxAF.R
#Max calculation

data <- reactive({
  myPrev = 1/input$prev
  if(input$inh=="monoallelic"){
    myMaxAF = (1/2) * myPrev * input$hetA * input$hetG * (1/input$pen)
  }
  if(input$inh=="biallelic"){
    myMaxAF = sqrt(myPrev) * input$hetA * sqrt(input$hetG) * (1/sqrt(input$pen))
  }
  myMaxAC = qpois(p=as.numeric(input$CI),
                  lambda=(input$popSize)*(myMaxAF))
  return(list(myMaxAF,myMaxAC))
})

#output$maxAF <- renderText({signif(data()[[1]],3)})

uncover_maxaf <- reactive ({

  condformat(intBED()) %>%
    rule_fill_discrete(ClinVar, expression= ClinVar !=".",
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(MutationAssessor, expression= MutationAssessor =='H',
                       colours= c("TRUE"="red", "FALSE"="green")) %>%
    rule_fill_discrete(c(M_CAP, AF_gnomAD),
                       expression= M_CAP > 0.025 & AF_gnomAD <0.05,
                       colours= c("TRUE"= "red", "FALSE"= "green")) %>%
    rule_fill_discrete(c(start, end),
                       expression= MutationAssessor =='H' &
                         ClinVar !="." & M_CAP != "TRUE" &
                         AF_gnomAD < (signif(data()[[1]],3)) ,
                       colours= c("TRUE"= "yellow", "FALSE"= ""))%>%
    rule_css(c(start, end),
             expression= ifelse(grepl("H|M", MutationAssessor) &
                                  ClinVar !="." & M_CAP != "TRUE" &
                                  AF_gnomAD < (signif(data()[[1]],3)),
                                "red", "green"), #& CADD_PHED > 20
             css_field = "color")

})


#server-binomial.R
#####output binomial distribution#####
df_subset <- reactive({
  validate(
    need(input$start_gp != "" , "Please, select genomic position"))
  mysample()
  validate(need(
    input$start_gp %in% mysample()$start,
    "Please, select a valide genomic position"))

  a <- subset(mysample(),start == as.numeric(input$start_gp) &
                end == as.numeric(input$start_gp))
  cov_ex<- as.numeric(a$coverage)
  print(cov_ex)
})


#ui.R
#' @title unCOVERApp user interface
#'
#' @description This function allows you to open graphical interface.
#' @param agree visualization.
#' @keywords app
#' @export
#' @examples
#' uncoverappLib.run()

require(shiny)
require(shinyWidgets)
require(shinyBS)
require(shinyjs)
require(markdown)
require(DT)
#options(repos = BiocInstaller::biocinstallRepos())
#getOption("repos")
#require(data.table)
#require(dplyr)
require(Gviz)
require(Homo.sapiens)
require(OrganismDbi)
require(stringr)
require(condformat)
require(shinyjs)
require(shinycssloaders)
require(bedr)
require(Rsamtools)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
#require(BSgenome.Hsapiens.UCSC.hg19)
#require(BSgenome.Hsapiens.UCSC.hg38)
require(EnsDb.Hsapiens.v75)
require(EnsDb.Hsapiens.v86)
require(org.Hs.eg.db)
options("browser" = "xdg-open")

first.file <- system.file(
  "extdata",
  "intro.md",
  package = "uncoverappLib"
)
second.file <- system.file(
  "extdata",
  "script.js",
  package = "uncoverappLib"
)
third.file <- system.file(
  "extdata",
  "CONTACTS.md",
  package = "uncoverappLib"
)



intro_processing <- system.file(
  "extdata",
  "prep_input.md",
  package = "uncoverappLib"
)

mdStyle <- "margin-left: 30px; margin-right: 30px"


intro <- function() {
  tabPanel("Home",
           shiny::includeMarkdown(first.file),
           #includeMarkdown(file.path(".","intro.md")),
           style=mdStyle
  )
}



preprocess <- function() {
  shiny::tabPanel("Processing and Statistical Summary",
                   shiny::includeMarkdown(intro_processing),
                       #shiny::includeMarkdown(file.path(".", "prep_input.md")),
                       h1(strong("Prepare your input file")),
                       fluidPage(
                         sidebarLayout(
                      sidebarPanel(
                        shiny::selectInput("Genome",
                                           label = "Reference Genome",
                                           choices = c("hg19",
                                                       "hg38"),
                                           selected = "UCSC genome"),
                        hr(),
                        shinyWidgets::pickerInput("notation",
                                                  label = "Chromosome Notation",
                                                  choices = c("chr", "number"),
                                                  options = list(`actions-box` = TRUE),
                                                  multiple =FALSE),

                        shinyWidgets::pickerInput("MAPQ",
                                                  label = "Minum Mapping Quality (MAPQ)",
                                                  choices = c(1:1000),
                                                  options = list(`actions-box` = TRUE), multiple =FALSE),

                        shinyWidgets::pickerInput("base_qual",
                                                  label = "Minimum Base Quality",
                                                  choices = c(1:1000),
                                                  options = list(`actions-box` = TRUE), multiple =FALSE),


                        # sliderInput("min_coverage",
                        #             label = "Minimum Coverage Depth (BED only)",
                        #             value = 20,
                        #             min = 1,
                        #             max = 1000,
                        #             step = 1)
                        hr(),

                        fileInput(inputId = "gene1",
                                  label = "Load input filtering file",
                                  accept = c("text/csv",
                                             ".zip",
                                             ".gz",
                                             ".bed",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        helpText("Choose the type of your input file"),
                        radioButtons("type_file", "Type of file",
                                     choices = c("List of genes name" = "gene_name",
                                                 "Target Bed "= "target_bed"),
                                     selected = "gene_name"),
                        br(),
                        br(),
                        helpText("The file must contain one path per row"),
                        fileInput(inputId = "bam_list",
                                  label = "Load your genomic data(BAM list or BED coverage list)",
                                  accept = c("text/csv",
                                             ".zip",
                                             ".gz",
                                             ".bed",
                                             "text/comma-separated-values,text/plain",
                                             ".list")),
                        helpText("Choose the type of your input file"),
                        radioButtons("type_coverage", "File format",
                                     choices = c("BAM file" = "bam",
                                                 "BED coverage file"= "bed"),
                                     selected = "bed"),
                        hr(),
                        hr(),
                        downloadButton("summary", "Statistical_Summary", class = "btn-primary",
                                       style='color: #fff; background-color: #27ae60;
                                       border-color: #fff;
                                       padding: 15px 14px 15px 14px;
                                       margin: 15px 5px 5px 5px; ')
                      ),
                      shiny::mainPanel(
                        tabsetPanel(
                          tabPanel(title= "input for uncoverapp",
                                   shinycssloaders::withSpinner(
                                     DT::DTOutput()("input1"))
                          )
                        ) )
                    )))
}


myHome <- function() {
  tabPanel("Coverage Analysis",
           h1(strong("Interactive web-application to visualize and
                     annotate low-coverage positions in clinical sequencing")),
           #titlePanel("Coverage sequencing Explorer"),
           helpText(em("Note:Select input options",
                       span("Upload your input bed
                            file with columns: chromosome, start, end, coverage
                            by sample and nucleotide count", style = "color:blue"))),
           shinyjs::useShinyjs(),
           includeScript(second.file),
           sidebarPanel(
             selectInput("UCSC_Genome",
                         label = "Reference Genome",
                         choices = c("hg19",
                                     "hg38"),
                         selected = "UCSC genome"),
             hr(),

             textInput(inputId = "Gene_name",
                       label= "Gene name"),
             #actionButton("button1",label= "Apply"),
             shinyBS::bsButton("button1",label= "Apply",  icon = icon("power-off"),
                               style = "success", size = "extra-small"),
             shinyjs::hidden(p(id = "text1", "Running.....")),
             #actionButton("remove",label= "Refresh"),
             shinyBS::bsButton("remove",label= "Refresh",  icon = icon("power-off"),
                               style = "success", size = "extra-small"),
             helpText(em("write gene name and push apply button")),
             hr(),

             textInput(inputId ="Chromosome",
                        label = "Chromosome"),

             #shinyWidgets::pickerInput("Chromosome",
              #                         label = "Chromosome",
               #                        choices = c("chr1", "chr2","chr3", "chr4","chr5",
                #                                   "chr6","chr7","chr8","chr9","chr10",
                 #                                  "chr11","chr12","chr13","chr14","chr15",
                  #                                 "chr16","chr17",
                   #                                "chr18","chr19", "chr20", "chr21",
                    #                               "chr22", "chrX", "chrY","chrM",
                     #                              names("file1")),
                      #                 options = list(`actions-box` = TRUE),
                       #                multiple =FALSE),
             hr(),
             shinyWidgets::pickerInput("coverage_co",
                                       label = "Coverage threshold",
                                       choices = c(1:1000, "all"),#names("file1"),
                                       options = list(`actions-box` = TRUE), multiple =FALSE),
             helpText(em("Select minimum value as coverage threshold")),
             hr(),
             textInput(inputId = "Sample",
                       label= "Sample"),

             helpText(em("Select name of sample for coverage analysis as reported
             in bam list input.
                         Example:example_POLG.bam")),

             hr(),
             splitLayout(cellWidths = c("30%", "70%"),
                         textInput(inputId = "transcript_id",
                                   label= "Transcript number"),
                         textInput(inputId = "id_t",
                                   label= "Transcript ID")),

             helpText(em("Retrieve your favourite transcript number from UCSC exons")),


             hr(),
             shinyWidgets::pickerInput("exon_number",
                                       label= "exon number",
                                       choices = c(1:150),
                                       options = list(`actions-box` = TRUE), multiple =FALSE),
             #actionButton("button5",label= "Make exon"),
             shinyBS::bsButton("button5",label= "Make exon",
                               icon = icon("power-off"), style = "success",
                               size = "extra-small"),
             shinyjs::hidden(p(id = "text1", "Running.....")),
             #actionButton("remove5",label= "Refresh"),
             shinyBS::bsButton("remove5",label= "Refresh",
                               icon = icon("power-off"), style = "default",
                               size = "extra-small"),
             helpText(em("zooming one exon")),
             hr(),
             hr(),

             textInput(inputId = "Start_genomicPosition",
                       label = "START genomic position"),

             textInput(inputId = "end_genomicPosition",
                       label = "END genomic position"),
             helpText(em("change genomic interval for zooming")),


             hr(),
             hr(),
             textInput(inputId = "query_Database",
                       label= "Region coordinates"),


             helpText(em("write to expand dbNSFP-annotated genomic positions.
                         For example 2:166845670-166930180")),
             hr(),

             hr(),
             downloadButton("downloadData", "Download", class = "btn-primary",
                            style='padding:4px; font-size:120%'),
             hr(),
             hr(),
             shiny::tags$button(
               id = 'close',
               type = "button",
               class = "btn action-button",
               style='color: white; background-color: #dd4b39;
               padding:4px; font-size:120%',
               onclick = "setTimeout(function(){window.close();},500);",
               # close browser
               "Close App")
            ),
           mainPanel(
             fileInput(inputId = "file1",
                       label = "Select input file",
                       accept = c("text/csv",
                                  ".bedGraph",
                                  ".bedGraph.gz",
                                  ".zip",
                                  ".gz",
                                  ".bed",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")),
             checkboxInput("header", "Header", TRUE),
             shinyBS::bsButton("pileup",
                               label= "load input file",
                               icon = icon("file"), style = "info",
                               size = "default"),
             #shinyBS::bsButton("example_data",
             #                 label= "load example dataset",
             #                icon = icon("file"), style = "info",
             #               size = "extra-small"),
             #helpText(em("load example dataset with
             #           base coverage counts of POLG gene")),
             hr(),

             tabsetPanel(
               tabPanel("bed file", shinycssloaders::withSpinner(
                 DT::DTOutput()("text"))),
               #tabPanel("UCSC gene",
               #         DT::DTOutput()("ccg")),
               #tabPanel("UCSC exons",
               #         helpText(em("Extract protein coding positions from
               #                     UCSC")), DT::DTOutput()("exon_pos")),
               tabPanel("Low-coverage positions",
                        DT::DTOutput()("text_cv")),
               #tabPanel("Gene coverage", shinycssloaders::withSpinner(
                 #plotOutput("all_gene")),
                 #DT::DTOutput()('df.l')),
               #tabPanel("Exon coverage", helpText(em("This is a Gviz function
               #and it plots exon with ideogram, genome coordinates,
               #coverage information, Ensembl and UCSC gene annotation.The
              #annotation for the databases are directly fetched from Ensemb and
                      #all tracks will be plotted in a 3' -> 5' direction.")),
                        #plotOutput("ens",
                                   #dblclick = "plot1_dblclick",
                                   #brush = brushOpts( id = "plot1_brush",
                                                      #resetOnNew = TRUE)),
                        #DT::DTOutput()("tabExon")),
               #tabPanel("Zoom to sequence",
                        #helpText(em("choose a few genomic intervals")),
                        #plotOutput("sequence")),
               tabPanel("Annotations on low-coverage positions",
                        helpText(em("dbSNP-annotation collects all
                                    consequences found in VEP-defined
                                    canonical transcripts")),
                        shinycssloaders::withSpinner(
                          condformat::condformatOutput("uncover_position"))),
               id = "tabSet"
             ),
             hr(),

             fluidRow(
               column(12,DT::DTOutput()("x4"))
             )
           )
  )
} 
myTab1 <- function() {
  tabPanel("Calculate AF by allele frequency app",

           # Application title
           titlePanel("Maximum credible population allele frequency"),

           ##### Bootstrap method for page costruction
           fluidPage(
             fluidRow(
               ##### Sidebar
               column(8,wellPanel(radioButtons("inh",
                                               "Inheritance:",
                                               choices = list("monoallelic",
                                                              "biallelic"),
                                               selected = "monoallelic"),

                                  numericInput("prev","Prevalence = 1 in ...
                                               (people)",
                                               min = 1,max = 1e8,value = 500),
                                  options = NULL),
                      br(),
                      sliderInput("hetA","Allelic heterogeneity:",
                                  min = 0, max = 1,value = 0.1),
                      sliderInput("hetG",
                                  "Genetic heterogeneity:",
                                  min = 0, max = 1,value = 1),
                      br(),
                      sliderInput("pen", "Penetrance:",
                                  min = 0, max = 1, value = 0.5))),
             br(),
             column(8,
                    h3("Maximum credible population AF:"),
                    h2(textOutput("maxAF"),align="center",style = "color:blue")),
             column(8,
                    h3("Uncover position",
                       helpText(em("Low-coverage positions excluding sites
                                   annotated as variants with AF> maxAF
                                   (default maxAF value: 5%)"),align="center",
                                style="color:blue"),
                       style = "font-size: 100%; width: 100%",
                       shinycssloaders::withSpinner(
                       condformat::condformatOutput("uncoverPosition")))),
             br(),
             br(),
             downloadButton("download_maxAF", "Download_maxAF",
                            class = "btn-primary",
                            style='padding:4px; font-size:80%',
                            helpText("download low coverage
                                     position dbSNFP-annotation filtered by
                                     maximum allele frequency",
                                     class = "btn-primary",
                                     style='padding:4px; font-size:60%'))
             #)
           ))}


myTab2 <- function() {
  tabPanel("Binomial distribution",
           titlePanel("Binomial distribution "),
           fluidRow(
             column(4,(numericInput("p",
                                    "Allele Fraction",
                                    min = 0,
                                    max = 1,
                                    value = 0.05)),
                    helpText(em("the expected fraction of variant reads
                    (probability of success)",
                    align="center",
                    style="color:gray")),
                    hr(),
                    numericInput("num_all",
                                 "Variant reads",
                                 min=0,
                                 max=100000,
                                 value=10),

                    helpText(em("the minimum number of variant reads required
                                by the user to support variant evidence,
                                (number of successes)"),align="center",
                             style="color:gray"),
                    hr(),

                    textInput(inputId = "start_gp",
                              label = "Genomic position"),

                    #textInput(inputId = "end_gp",
                     #         label = "END genomic position"),
                    helpText(em("Specify a genomic position of interest gene
                                  in which calculating the
                                binomial distribution"))),


             column(4,
                    h2("consideration:"),
                    h3(htmlOutput("ci"))),
             column(4,
                    h2("Binomial Distribution", plotOutput("bd"))),
             column(10, h2("Cumulative distribution function"),
                    h3(plotOutput("pbinom"))))

  )}


myabout <- function() {
  tabPanel("Contacts",
           includeMarkdown(third.file),
           style=mdStyle
  )
}

ui <- shinyUI(
  tagList(
    shiny::tags$head(
      shiny::tags$style(HTML("
        .navbar { background-color: #FFFFFF;height: 150px;}
        .navbar-brand img { height: 150px;margin: 0px; }
        .navbar-nav > li > a { color: #333; font-weight: bold; font-size: 15px; margin-top: 60px }
        .navbar-nav > li > a:hover { color: #007BFF; }
        .navbar-brand { padding: 0px 10px; }
      "))
    ),
    
    navbarPage(
      title = shiny::tags$a(
        href = "#",
        class = "navbar-brand",
        img(src = "logo.png")
      ),
      windowTitle = "uncoverApp",
      
      # Dropdown menu per alcune funzionalità principali
      navbarMenu("Menu",
                 intro(),
                 preprocess(),
                 myHome(),
                 myTab1(),
                 myTab2()
                 #myabout()
      ),
      
      # Tab separati
      #myTab1(),
      #myTab2(),
      myabout()
    )
  )
)

# Remove the extra closing brace that was causing the error
# The commented-out alternative UI definition can be removed or kept as a comment
#  ui <- shinyUI(navbarPage(
#    title = div(
#      img(src = "logo.png", 
#          height = "100px", 
#          style = "margin-top: 10px; margin-bottom: -10px;"),
#      style = "margin: 4; padding: 4;"
#    ),
#    windowTitle = "uncoverApp",  # This sets the browser tab title
#    intro(),
#    preprocess(),
#    myHome(),
#    myTab1(),
#    myTab2(),
#    myabout()
#  ))

#server.R
server <- function (input, output, session){
  options(shiny.maxRequestSize=30*1024^2)

  script1 <- system.file(
    "extdata",
    "Rpreprocessing.R",
    package = "uncoverappLib")


  #attach static scripts


  output$dependence = downloadHandler(filename="Rpreprocessing.R",
                                      content=function(file){
                                        file.copy(script1,file)
                                      })
  #source script to load dataset or example file
  source('server-preprocess.R', local= TRUE)
  output$input1<- renderDataTable({
    options(shiny.sanitize.errors = TRUE)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "table construction in progress",
                 detail = 'This may take a while', value = 0)
    Sys.sleep(0.1)
    validate(
      need(try(!is.null(coverage_input())), "please upload a file with HGNC
      gene names and absolute path(s) to BAM file"))
    coverage_input()
  })

  output$summary <- downloadHandler(
    filename = function() {
      paste('statistical_summary',input$file1,Sys.Date(), '.txt')
    },
    content = function(file){
      # Load randomizer data
      omim_gene <- read.table("/home/anna/Desktop/sys_ndd_2025_subset.tsv", 
                        header=TRUE, sep="\t", stringsAsFactors=FALSE)

      # Get only SYMBOL and MIM_Number columns  

      # Join with statistical summary (keep all from stat_summ)
      original_data <- stat_summ()
      joined_data <- merge(original_data, omim_gene, by="SYMBOL", all.x=TRUE)
      print(joined_data)

      # Write joined result
      write.table(joined_data, file, sep='\t', quote=F, row.names = F, col.names = T)
      }
    )

  source('server-reactiveDF.R', local= TRUE)

  output$text<- #DT::renderDataTable({
    DT::renderDataTable({
      #validate(need(input$file1 != "", "Please, upload your file"))
      validate(need(ncol(mydata()) != "0", "Please, upload your file"))
      mydata()})

  output$text_cv <- DT::renderDataTable({
    validate(need(input$Gene_name != "" & input$Sample !="",
                  "Please select all required input: Gene, Chromosome,
                  Coverage threshold and Sample"))
    filtered_low()})

  #source script to reactive gene and exon table
  source('server-tables.R', local= TRUE)

  #source script to plot all gene coverage obteined by gviz
  #plot and to table low coverage position in each exon
  source('server-plots.R', local=TRUE)

  output$all_gene<- renderPlot({
    validate(
      need(ncol(mydata()) != "0", "Unrecognized data set: Please
           upload your file"))  # check error message
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

  output$df.l<- DT::renderDataTable({
    table1()})

  #plot with sequence reference

  output$sequence<- renderPlot({
    validate(
      need(ncol(mydata()) != "0",
           "Unrecognized data set: Please load your file"))
    validate(
      need(input$Start_genomicPosition < input$end_genomicPosition,
           "Please selct the right genomic position: end position is
           lower than start "))
    options(shiny.sanitize.errors = TRUE)
    progress <- shiny::Progress$new()
    for (i in 1:40) {
      progress$set(message = "Please wait a few minutes: Making plot",
                   detail = 'This may take a while', value = i)
      Sys.sleep(1.0)}
    on.exit(progress$close())
    Sys.sleep(0.1)
    p3()
  })

  #source script to dbNSFP annotation
  source('server-annotation.R', local=TRUE)
  output$tabExon<- DT::renderDataTable({
    validate(
      need(input$Gene_name != "" & input$Sample != "",
           "Please select all required input: Gene, Chromosome,
           Coverage threshold and Sample")
    )
    validate(
      need(try(!is.null(intBED())),'Unrecognized coordinates:
      Please change exon number input and be sure that input box
      "Region coordinates" is filled.
           Please check your R environment has annotation file loaded.'))
    #print(try(is.null(intBED())))
    if (is.null(intBED()))
      return(NULL)
    a= intBED() %>%
      dplyr::select(MutationAssessor,  ClinVar)
    df=data.frame(length(which(a$ClinVar!='.')),
                  length(which(grepl("H|M", a$MutationAssessor))))
    colnames(df)= c('Low coverage positions in ClinVar',
                    'Low coverage positions with High or Medium impact')
    return(df)

  })

  #table output wiht Condformat pkg

  output$uncover_position<- condformat::renderCondformat({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "table construction in progress", value = 0)
    Sys.sleep(0.1)
    validate(
      need(input$Gene_name != "" & input$Sample != "",
           "Please select all required input: Gene, Chromosome,
           Coverage threshold and Sample"))
    validate(
      need(try(!is.null(intBED())),'Unrecognized coordinates:
      Please change exon number input and be sure that input box
      "Region coordinates" is filled. This function works for only Hg19 coordinate!
           Please check your R environment has annotation file loaded.')
    )
    if(is.null(condform_table()))
      return(NULL)
    #print(head(condform_table()))
    condform_table() %>%
     dplyr::select(seqnames, start, end, coverage, REF, ALT,
                   dbsnp, GENENAME, PROTEIN_ensembl,  MutationAssessor,SIFT,
                    Polyphen2,M_CAP,CADD_PHED,AF_gnomAD,
                    ClinVar,clinvar_MedGen_id,HGVSc_VEP,
                  HGVSp_VEP)
    })

  #download data wiht conditionalFormatting

  output$downloadData <- downloadHandler(
    filename = function() {
      paste('download',input$file1,Sys.Date(), '.xlsx', sep='')
    },
    content = function(file){
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Sheet1")
      negStyle <- openxlsx::createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
      posStyle <- openxlsx::createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
      highlighted<-openxlsx::createStyle (fgFill = "yellow")
      hs <- openxlsx::createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF",
                        fontSize=12,
                        fontName="Arial Narrow", fgFill = "#4F80BD")
      openxlsx::writeData(wb, "Sheet1",condform_table(), headerStyle = hs)
      openxlsx::conditionalFormatting(wb, "Sheet1",cols=19,
                            rows=(1:nrow(condform_table())+1),rule='$S2=="H"',
                            style =negStyle)
      openxlsx::conditionalFormatting(wb, "Sheet1",cols=19,
                            rows=(1:nrow(condform_table())+1),
                            rule='$S2!="H"', style =posStyle)
      openxlsx::conditionalFormatting(wb, "Sheet1",cols=22,
                            rows=(1:nrow(condform_table())+1),
                            rule='$V2=="D"', style =negStyle)
      openxlsx::conditionalFormatting(wb, "Sheet1",cols=22,
                            rows=(1:nrow(condform_table())+1),
                            rule='$V2!="D"', style =posStyle)
      openxlsx::conditionalFormatting(wb, "Sheet1", cols=23,
                            rows=(1:nrow(condform_table())+1),
                            rule=">=20", style = negStyle)
      openxlsx::conditionalFormatting(wb, "Sheet1", cols=23,
                            rows=(1:nrow(condform_table())+1),
                            rule="<20", style = posStyle)
      openxlsx::conditionalFormatting(wb, "Sheet1", cols=24,
                            rows=(1:nrow(condform_table())+1),
                            rule="<=0.05", style = negStyle)
      openxlsx::conditionalFormatting(wb, "Sheet1", cols=25,
                            rows=(1:nrow(condform_table())+1),
                            rule='!="."', style = negStyle)
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )


  # source code for calculator maxAF

  source('server-maxAF.R', local=TRUE)
  output$maxAF <- renderText({signif(data()[[1]],3)})
  output$uncoverPosition <- condformat::renderCondformat ({
    validate(
      need(ncol(mydata()) != "0",
           "Unrecognized data set: Please load your file"),
      need(input$Gene_name != "" &
             input$Sample != "", "Please select all required input: Gene,
           Chromosome, Coverage threshold and Sample")
    )
    validate(
      need(try(!is.null(uncover_maxaf())),
           'Unrecognized coordinates: Please change exon number input and
           be sure that input box "Region coordinates" is filled')
    )
    uncover_maxaf()
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "table construction in progress", value = 0)
    Sys.sleep(0.1)

    uncover_maxaf() %>%
      dplyr::select(seqnames, start, end, coverage, counts,REF,
                    ALT, dbsnp, GENENAME, PROTEIN_ensembl, MutationAssessor,
                    SIFT,Polyphen2,M_CAP,CADD_PHED,AF_gnomAD,
                    ClinVar,clinvar_MedGen_id,
                    HGVSc_VEP,HGVSp_VEP)
  })

  ###download excel with maxAF###

  output$download_maxAF <- downloadHandler(
    filename = function() {
      paste('download_uncover_maxAF',Sys.Date(), '.xlsx', sep='')
    },
    content = function(file){
      wb1 <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb1, "Sheet2")
      negStyle <- openxlsx::createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
      posStyle <- openxlsx::createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
      highlighted<-openxlsx::createStyle (fgFill = "yellow")
      hs <- openxlsx::createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF",
                        fontSize=12,
                        fontName="Arial Narrow", fgFill = "#4F80BD")
      openxlsx::writeData(wb1, "Sheet2",uncover_maxaf(), headerStyle = hs)
      openxlsx::conditionalFormatting(wb1, "Sheet2",cols=19,
                            rows=(1:nrow(uncover_maxaf())+1),
                            rule='$S2=="H"', style =negStyle)
      openxlsx::conditionalFormatting(wb1, "Sheet2",cols=19,
                            rows=(1:nrow(uncover_maxaf())+1),
                            rule='$S2!="H"', style =posStyle)
      openxlsx::conditionalFormatting(wb1, "Sheet2",cols=22,
                            rows=(1:nrow(uncover_maxaf())+1),
                            rule='$V2=="D"', style =negStyle)
      openxlsx::conditionalFormatting(wb1, "Sheet2",cols=22,
                            rows=(1:nrow(uncover_maxaf())+1),
                            rule='$V2!="D"', style =posStyle)
      openxlsx::conditionalFormatting(wb1, "Sheet2", cols=23,
                            rows=(1:nrow(uncover_maxaf())+1),
                            rule=">=20", style = negStyle)
      openxlsx::conditionalFormatting(wb1, "Sheet2", cols=23,
                            rows=(1:nrow(uncover_maxaf())+1),
                            rule="<20", style = posStyle)
      openxlsx::conditionalFormatting(wb1, "Sheet2", cols=24,
                            rows=(1:nrow(uncover_maxaf())+1),
                            rule="<=0.05", style = negStyle)
      openxlsx::conditionalFormatting(wb1, "Sheet2", cols=25,
                            rows=(1:nrow(uncover_maxaf())+1),
                            rule='!="."', style = negStyle)
      openxlsx::saveWorkbook(wb1, file, overwrite = TRUE)
      #condformat2excel(condform_table(), file, sheet_name = "Sheet1",
      #                overwrite_wb = F, overwrite_sheet = T)
    }
  )

  #source script for binomial distribution given genomic position
  source('server-binomial.R', local=TRUE)
  output$bd<- renderPlot({
    library(stats)
    library(stats)
    p=seq(0,1,by=0.0001)

    ic<-qbinom(p = c(0.025, 0.975), size = df_subset(), prob =input$p)
    barplot(dbinom(x=1:df_subset(), size=df_subset(), p= input$p) ,
            main = paste(ic))
  })

  output$pbinom<- renderPlot({
    library(stats)
    #c=input$num_all
    p=seq(0,1,by=0.0001)
    Fx=pbinom(q=1:df_subset(), size=df_subset(), prob = input$p)
    npro=Fx[as.numeric(input$num_all)]
    leg= (1-npro)*100

    #leggend= (1-Fx[input$num_all]) *100
    #print(leggend)
    if (is.na(npro)){
      txt= "under threshold"}
    else{
      txt= paste("the probability of detecting more of",
                 input$num_all, "reads wiht variant alleles would be:",
                 leg, "%")
    }
    plot(1:df_subset(), Fx, type='h',xlab= "number of trials",
         main= txt)
    abline(v=input$num_all, col= "red")
  })

  output$ci<- renderText({
    #return(df_subset())
    ci_print<-qbinom(p = c(0.025, 0.975), size = df_subset(), prob =input$p)
    ci_in=c(ci_print[1]:ci_print[2])
    number= ci_print[1]
    number2=ci_print[2]
    thr= as.numeric(df_subset())
    # if (thr < input$coverage_co) {
    if (number2 < input$num_all) {
      print(paste("<span style=\"color:red\">according to the binomial probability model,
                  there is 95% probability to observe from </span>", number,
                  "<span style=\"color:red\">to</span>",
                  number2,
                  "<span style=\"color:red\"> variant-supporting reads .
                  </span>"))

    }else{
      print(paste("<span style=\"color:blue\">according to the binomial probability model,
                  there is 95% probability to observe from </span>", number,
                  "<span style=\"color:blue\">to</span>",
                  number2,
                  "<span style=\"color:blue\"> variant-supporting reads
                  .</span>"))
    }
  })
  observe({
    if (input$close > 0) stopApp()
  })
}

#server-reactiveDF.R
###make reactive database given reference genome
txdb= reactive({if (input$UCSC_Genome == "hg19"){
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene}
  else{
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene}
})

#################make reactive dataframe given input file

mydata <- reactiveVal()
observeEvent(input$file1, {

  # Read uploaded table (user says this is already like a data frame)
  tmp <- read.table(input$file1$datapath,
                    header = input$header, stringsAsFactors = FALSE)

  # Ensure first three columns are named consistently
  if (ncol(tmp) >= 3) {
    colnames(tmp)[1:3] <- c("chromosome","start","end")
    head(tmp)
  }

  # If the uploaded table already has sample_* columns, just use it as-is
  if (any(grepl("^sample_", names(tmp)))) {
    mydata(tmp)
    return(invisible(NULL))
  }

  # Otherwise, if it has a 4th column as coverage, map it to one sample
  sample_label <- tools::file_path_sans_ext(basename(input$file1$name))
  cov <- if (ncol(tmp) >= 4) suppressWarnings(as.integer(tmp[[4]])) else rep(NA_integer_, nrow(tmp))
  new_df <- data.frame(
    chromosome = tmp$chromosome,
    start = as.integer(tmp$start),
    end = as.integer(tmp$end),
    cov = cov,
    stringsAsFactors = FALSE
  )
  colnames(new_df)[4] <- paste0("sample_", sample_label)

  mydata(new_df)
})

observeEvent(input$pileup, {

  tmp_pileup <- coverage_input()
  colnames(tmp_pileup)[1:3] <- c("chromosome","start","end")

  # Determine how many coverage columns we have and how many samples
  data_cols <- setdiff(colnames(tmp_pileup), c("chromosome","start","end"))
  ncols <- length(data_cols)
  samples <- name_sample()
  k <- length(samples)

  if (ncols == k && k > 0) {
    # One coverage column per sample: name them as sample_* only
    colnames(tmp_pileup)[-1:-3] <- paste0("sample_", samples)
  } else if (ncols == 2 * k && k > 0) {
    # Two columns per sample (legacy BAM): interleave as sample_* and nucleotide_*
    interleaved <- as.vector(rbind(paste0("sample_", samples), paste0("nucleotide_", samples)))
    colnames(tmp_pileup)[-1:-3] <- interleaved
  } else {
    # Fallback: generic sample_1..sample_n
    colnames(tmp_pileup)[-1:-3] <- paste0("sample_", seq_len(ncols))
  }

  # Enforce types; chromosome already normalized upstream in coverage_input()
  tmp_pileup$chromosome <- as.character(tmp_pileup$chromosome)
  tmp_pileup$start <- as.integer(tmp_pileup$start)
  tmp_pileup$end <- as.integer(tmp_pileup$end)

  mydata(tmp_pileup)
})

mysample<-reactive({
  dat <- mydata()
  if (is.null(dat)) return(NULL)
  sel <- input$Sample
  if (is.null(sel) || sel == "") return(NULL)

  # Accept either raw sample id, with/without extension, or full column name
  candidates <- unique(c(
    sel,
    sub("^sample_", "", sel),
    tools::file_path_sans_ext(sel),
    paste0("sample_", sel),
    paste0("sample_", tools::file_path_sans_ext(sel))
  ))
  col <- intersect(paste0("sample_", candidates), names(dat))
  if (length(col) == 0) return(NULL)
  col <- col[1]

  out <- dat[, c("chromosome","start","end", col), drop = FALSE]
  names(out)[4] <- "coverage"
  out
})

filtered_low<- reactive ({

  if (is.null(mysample()))
    return(NULL)
  df <- mysample()
  thr <- input$coverage_co
  if (identical(thr, "all")) {
    dplyr::filter(df, chromosome == Chromosome())
  } else {
    dplyr::filter(df, chromosome == Chromosome(), coverage <= as.numeric(thr))
  }

})

filtered_high<- reactive ({
  if (is.null(mysample()))
    return(NULL)
  df <- mysample()
  thr <- input$coverage_co
  if (identical(thr, "all")) {
    dplyr::filter(df, chromosome == Chromosome())
  } else {
    dplyr::filter(df, chromosome == Chromosome(),
                  coverage > as.numeric(thr))
  }
})
