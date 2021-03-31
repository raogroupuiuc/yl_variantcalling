' Analysis of actual data 
Started March 4, 2021 from workshop practised code 
used for Xylose mutants 
Step 2- analyse the imported vcf files '

getwd()
setwd("D:/Box Sync/02 - RNA seq work/2101-YL-Evolution (Sangdo)/Git_shared/yl_variantcalling/R_scripts/")

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

rm(list = ls())

load("01_dataimport.Rdata")

head(rowRanges(readvcf.plasmid))

head(info(readvcf.plasmid))
colnames(info(readvcf.plasmid))
# [1] "AC"              "AF"              "AN"              "BaseQRankSum"    "DP"              "DS"              "END"            
# [8] "ExcessHet"       "FS"              "InbreedingCoeff" "MLEAC"           "MLEAF"           "MQ"              "MQRankSum"      
# [15] "QD"              "RAW_MQandDP"     "ReadPosRankSum"  "SOR"  

info(readvcf.plasmid)$DP # use this to make histograms 

geno(readvcf.plasmid)
# List of length 8
# names(8): GT AD DP GQ MIN_DP PL RGQ SB

head(geno(readvcf.plasmid)$DP)

head(rowRanges(readvcf.strains))

head(info(readvcf.strains))
colnames(info(readvcf.strains))
# [1] "AC"              "AF"              "AN"              "BaseQRankSum"    "DP"              "DS"              "END"            
# [8] "ExcessHet"       "FS"              "InbreedingCoeff" "MLEAC"           "MLEAF"           "MQ"              "MQRankSum"      
# [15] "QD"              "RAW_MQandDP"     "ReadPosRankSum"  "SOR"  

info(readvcf.strains)$DP # use this to make histograms 

geno(readvcf.strains)
# List of length 8
# names(8): GT AD DP GQ MIN_DP PL RGQ SB

head(geno(readvcf.strains)$DP)

# write.table(geno(readvcf.plasmid)$DP, "results/readdepth-plasmid.txt", sep = "\t")
# write.table(geno(readvcf.strains)$DP, "results/readdepth-strains.txt", sep = "\t")
#  
# write.table(geno(readvcf.plasmid)$GT, "results/mutations-plasmid.txt", sep = "\t")
# write.table(geno(readvcf.strains)$GT, "results/mutations-strains.txt", sep = "\t")


# Working with genome annotations #### 

sig.hits.plasmid <- rowRanges(readvcf.plasmid)
sig.hits.strains <- rowRanges(readvcf.strains)

my_ld <- 1000

search.plasmid <- flank(sig.hits.plasmid, my_ld, both = TRUE)
search.strains <- flank(sig.hits.strains, my_ld, both = TRUE)

hits.plasmid <- findOverlaps(search.plasmid, gtf.gene.plasmid)

hits.plasmid2 <- data.frame(SNP = names(search.plasmid)[queryHits(hits.plasmid)],
                       Gene = gtf.gene.plasmid$ID[subjectHits(hits.plasmid)])

hits.plasmid3 <- findOverlaps(sig.hits.plasmid, gtf.gene.plasmid)

hits.plasmid4 <- data.frame(SNP = names(search.plasmid)[queryHits(hits.plasmid3)],
                            Gene = gtf.gene.plasmid$ID[subjectHits(hits.plasmid3)])
hits.plasmid4

hits.strains <- findOverlaps(search.strains, gtf.gene.strains)

hits.strains2 <- data.frame(SNP = names(search.strains)[queryHits(hits.strains)],
                            Gene = gtf.gene.strains$ID[subjectHits(hits.strains)])

hits.strains3 <- findOverlaps(sig.hits.strains, gtf.gene.strains)

hits.strains4 <- data.frame(SNP = names(search.strains)[queryHits(hits.strains3)],
                            Gene = gtf.gene.strains$ID[subjectHits(hits.strains3)])
hits.strains4
# 
# write.table(hits.plasmid4, "results/snpsingene-annotated-plasmid.txt", sep = "\t")
# write.table(hits.strains4, "results/snpsingene-annotated-strains.txt", sep = "\t")
# 

save.image("02-analysis-part1.Rdata")

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