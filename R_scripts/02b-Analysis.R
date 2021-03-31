' Analysis of actual data 
Started March 4, 2021 from workshop practised code 
used for Xylose mutants 
Step 2- analyse the imported vcf files '

getwd()
setwd("D:/Box Sync/02 - RNA seq work/2101-YL-Evolution (Sangdo)/Git_shared/yl_variantcalling/R_scripts//")

library(GenomicRanges)
library(GenomicFeatures)
library(Biostrings)
library(Rsamtools)
library(SummarizedExperiment)
library(magrittr)
library(VariantAnnotation)
library(snpStats)
library(ggplot2)
library(scater)


load("02-analysis-part1.Rdata")

# Identifying functional consequences of SNPs #### 

snps2search.plasmid <- subsetByOverlaps(sig.hits.plasmid, gtf.gene.plasmid)
snps2search.strains <- subsetByOverlaps(sig.hits.strains, gtf.gene.strains)

# genome.plasmid <- FaFile("data/pINT03-XYL123.fa")
txdb.strains <- makeTxDbFromGFF("data/GCA_009372015.1_YarliW29_genomic.gff",
                           format = "gff", organism = "Yarrowia lipolytica",
                           dataSource = "NCBI") # works

# Parents missing; import as GRanges and modify
pre.txdb.strains <- rtracklayer::import("data/GCA_009372015.1_YarliW29_genomic.gff")
pre.txdb.plasmid <- rtracklayer::import("data/pINT03-XYL123.gff")

# First make IDs unique
pre.txdb.plasmid$ID_orig <- pre.txdb.plasmid$ID
pre.txdb.plasmid$ID <- paste(pre.txdb.plasmid$type, pre.txdb.plasmid$ID_orig, sep = "-")

# Check how parents are coded in the one that works
table(pre.txdb.strains$type)
table(pre.txdb.plasmid$type)
pre.txdb.strains[pre.txdb.strains$type == "gene"] # no parent
pre.txdb.strains[pre.txdb.strains$type == "mRNA"] # parent is gene
pre.txdb.strains[pre.txdb.strains$type == "exon"] # parent is mRNA
pre.txdb.strains[pre.txdb.strains$type == "CDS"] # parent is mRNA

# Code parents for plasmid features
plasmid.parents <- rep(NA_character_, length(pre.txdb.plasmid))
plasmid.parents[pre.txdb.plasmid$type == "mRNA"] <- paste0("gene-", pre.txdb.plasmid$ID_orig[pre.txdb.plasmid$type == "mRNA"])
plasmid.parents[pre.txdb.plasmid$type == "exon"] <- paste0("mRNA-", pre.txdb.plasmid$ID_orig[pre.txdb.plasmid$type == "exon"])
plasmid.parents[pre.txdb.plasmid$type == "CDS"] <- paste0("mRNA-", pre.txdb.plasmid$ID_orig[pre.txdb.plasmid$type == "CDS"])

plasmid.parents <- as(plasmid.parents, "CharacterList")

for(i in seq_along(plasmid.parents)){
  if(all(is.na(plasmid.parents[[i]]))){
    plasmid.parents[[i]] <- character(0)
  }
}

pre.txdb.plasmid$Parent <- plasmid.parents

# Make TxDb from GRanges

txdb.plasmid <- makeTxDbFromGRanges(pre.txdb.plasmid)
# Warning message:
#   In .get_cds_IDX(mcols0$type, mcols0$phase) :
#   The "phase" metadata column contains non-NA values for features of type gene, mRNA, exon.
# This information was ignored.


seqinfo(snps2search.strains)
# Seqinfo object with 168 sequences from an unspecified genome:
seqinfo(gtf.gene.strains)
# Seqinfo object with 168 sequences from an unspecified genome; no seqlengths:
seqinfo(txdb.strains)
# Seqinfo object with 156 sequences from an unspecified genome; no seqlengths:
seqinfo(genome.strains)
# Seqinfo object with 168 sequences from an unspecified genome:

seqinfo(snps2search.plasmid)
# Seqinfo object with 168 sequences from an unspecified genome:
seqinfo(gtf.gene.plasmid)
# Seqinfo object with 169 sequences from an unspecified genome; no seqlengths:
seqinfo(txdb.plasmid)
# Seqinfo object with 1 sequence from an unspecified genome; no seqlengths:
seqinfo(genome.plasmid)
# Seqinfo object with 169 sequences from an unspecified genome:

# Hello errors #### 

rowRanges(snps2search.plasmid)
# Error in MatrixGenerics:::.load_next_suggested_package_to_search(x) : 
#   Failed to find a rowRanges() method for GRanges objects.

rowRanges(snps2search.strains)
# Error in MatrixGenerics:::.load_next_suggested_package_to_search(x) : 
#   Failed to find a rowRanges() method for GRanges objects.

pc.strains <- predictCoding(snps2search.strains, txdb.strains, genome.strains)
# Error in (function (classes, fdef, mtable)  : 
#             unable to find an inherited method for function 'predictCoding' for signature '"GRanges", "TxDb", "FaFile", "missing"'

pc.plasmid <- predictCoding(snps2search.plasmid, txdb.plasmid, genome.plasmid)

# table(pc$CONSEQUENCE)
# my_txnames <- select(my_txdb, pc$TXID, c("TXID", "TXNAME"), "TXID")
# head(my_txnames)
# identical(my_txnames$TXID, as.integer(pc$TXID))
# pc$TXNAME <- my_txnames$TXNAME
# 
# df <- data.frame(SNP = names(pc),
#                  Chromosome = seqnames(pc),
#                  Position = start(pc),
#                  Transcript = pc$TXNAME,
#                  Consequence = pc$CONSEQUENCE,
#                  AA_position = unlist(pc$PROTEINLOC),
#                  Ref_AA = pc$REFAA,
#                  Var_AA = pc$VARAA)
# write.csv(df, file = "protein_coding_variants.csv", row.names = FALSE)



# Converting from VCF to snpMatrix #### 

# Make a SNP matrix the DIY way if SnpMatrix is being a pain 
# (maybe it doesn't work on haploid data)

gt.strains <- geno(readvcf.strains)$GT
gt.plasmid <- geno(readvcf.plasmid)$GT

dim(gt.strains)
# 143 5
dim(gt.plasmid)
# 171 3 

table(gt.strains)
#  .   0   1   2   3   4 
# 18  90 558  41   2   1
table(gt.plasmid)
# .   0   1   2 
# 15  55 424  16 

colnames(gt.strains)
colnames(gt.plasmid)

mean(apply(gt.strains, 1, function(x) any(x %in% c("2", "3", "4")))) 
# 16.1% of loci have multiple alleles
mean(apply(gt.plasmid, 1, function(x) any(x %in% c("2")))) 
# 7.6% of loci have multiple alleles

# Function to make a binary matrix indicating presence/absence of each allele in the dataset
haploidSnpMat <- function(vcf){
  gt <- geno(vcf)$GT
  nalt <- lengths(rowRanges(vcf)$ALT)
  al.ind <- rep(seq_along(rowRanges(vcf)), times = nalt)
  al.names <- paste0(gsub("/.*$", "/", names(rowRanges(vcf)))[al.ind],
                     unlist(rowRanges(vcf)$ALT))
  out <- matrix(NA_integer_, nrow = length(al.names), ncol = ncol(gt),
                dimnames = list(al.names, colnames(gt)))
  curr.row <- 1L
  for(i in seq_len(nrow(gt))){
    theseals <- unique(gt[i,])
    theseals <- theseals[!is.na(theseals)]
    theseals <- theseals[theseals != "."]
    for(j in seq_len(nalt[i])){
      out[curr.row, which(gt[i,] %in% setdiff(theseals, as.character(j)))] <- 0L
      out[curr.row, which(gt[i,] == as.character(j))] <- 1L
      curr.row <- curr.row + 1L
    }
  }
  return(out)
}

mat.strains <- haploidSnpMat(readvcf.strains)
mat.plasmid <- haploidSnpMat(readvcf.plasmid)

dim(mat.strains)
# 170*5 = 850
dim(mat.plasmid)
# 185*3 = 555

is.na(mat.strains) %>% sum()
# 33
is.na(mat.plasmid) %>% sum()
#25 

table(mat.strains)
# 0   1 
# 215 602 = 817
table(mat.plasmid)
# 0   1 
# 90 440 = 530

# Basic statistics and filtering on SNPs #### 

summary(mat.plasmid)
summary(mat.strains)

save.image("02-analysis-part2.Rdata")

write.table(mat.plasmid, "results/snps-plasmid.txt", sep = "\t")
write.table(mat.strains, "results/snps-strains.txt", sep = "\t")

sessionInfo()
# R version 4.0.5 (2021-03-31)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
# [3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] scater_1.18.6               SingleCellExperiment_1.12.0 ggplot2_3.3.3              
# [4] snpStats_1.40.0             Matrix_1.3-2                survival_3.2-10            
# [7] VariantAnnotation_1.36.0    magrittr_2.0.1              SummarizedExperiment_1.20.0
# [10] MatrixGenerics_1.2.1        matrixStats_0.58.0          Rsamtools_2.6.0            
# [13] Biostrings_2.58.0           XVector_0.30.0              GenomicFeatures_1.42.2     
# [16] AnnotationDbi_1.52.0        Biobase_2.50.0              GenomicRanges_1.42.0       
# [19] GenomeInfoDb_1.26.4         IRanges_2.24.1              S4Vectors_0.28.1           
# [22] BiocGenerics_0.36.0        
# 
# loaded via a namespace (and not attached):
#   [1] viridis_0.5.1             httr_1.4.2                BiocSingular_1.6.0       
# [4] viridisLite_0.3.0         bit64_4.0.5               splines_4.0.5            
# [7] DelayedMatrixStats_1.12.3 scuttle_1.0.4             assertthat_0.2.1         
# [10] askpass_1.1               BiocFileCache_1.14.0      blob_1.2.1               
# [13] vipor_0.4.5               BSgenome_1.58.0           GenomeInfoDbData_1.2.4   
# [16] progress_1.2.2            pillar_1.5.1              RSQLite_2.2.5            
# [19] lattice_0.20-41           beachmat_2.6.4            glue_1.4.2               
# [22] colorspace_2.0-0          XML_3.99-0.6              pkgconfig_2.0.3          
# [25] biomaRt_2.46.3            zlibbioc_1.36.0           purrr_0.3.4              
# [28] scales_1.1.1              BiocParallel_1.24.1       tibble_3.1.0             
# [31] openssl_1.4.3             generics_0.1.0            ellipsis_0.3.1           
# [34] cachem_1.0.4              withr_2.4.1               crayon_1.4.1             
# [37] memoise_2.0.0             fansi_0.4.2               xml2_1.3.2               
# [40] beeswarm_0.3.1            tools_4.0.5               prettyunits_1.1.1        
# [43] hms_1.0.0                 lifecycle_1.0.0           stringr_1.4.0            
# [46] munsell_0.5.0             irlba_2.3.3               DelayedArray_0.16.3      
# [49] compiler_4.0.5            rsvd_1.0.3                rlang_0.4.10             
# [52] grid_4.0.5                RCurl_1.98-1.3            BiocNeighbors_1.8.2      
# [55] rappdirs_0.3.3            bitops_1.0-6              gtable_0.3.0             
# [58] DBI_1.1.1                 curl_4.3                  R6_2.5.0                 
# [61] gridExtra_2.3             GenomicAlignments_1.26.0  dplyr_1.0.5              
# [64] rtracklayer_1.49.5        fastmap_1.1.0             bit_4.0.4                
# [67] utf8_1.2.1                ggbeeswarm_0.6.0          stringi_1.5.3            
# [70] Rcpp_1.0.6                vctrs_0.3.6               sparseMatrixStats_1.2.1  
# [73] dbplyr_2.1.0              tidyselect_1.1.0