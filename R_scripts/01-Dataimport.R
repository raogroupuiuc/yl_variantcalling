' Analysis of actual data 
Started March 4, 2021 from workshop practised code 
used for Xylose mutants '

getwd()
setwd("D:/Box Sync/02 - RNA seq work/2101-YL-Evolution (Sangdo)/R analysis/")

# setwd("//DESKTOP-740VTGD/Desktop_work/Box Sync/02 - RNA seq work/2101-YL-Evolution (Sangdo)/R analysis")

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

# Reference genome #### 

ref.strains <- "data/GCA_009372015.1_YarliW29_genomic.fna"
# indexFa(ref.strains)
genome.strains <- FaFile(ref.strains)
seqinfo(genome.strains)
seqnames(seqinfo(genome.strains))

ref.plasmid <- "data/PO1f_plasmid.fna"
# indexFa(ref.plasmid)
genome.plasmid <- FaFile(ref.plasmid)
seqinfo(genome.plasmid)
seqnames(seqinfo(genome.plasmid))

# Import gene annotations #### 

gtf0.strains <- rtracklayer::import("data/GCA_009372015.1_YarliW29_genomic.gff")
gtf0.plasmid <- rtracklayer::import("data/PO1f_plasmid.gff")

gtf.gene.strains <- gtf0.strains[gtf0.strains$type == "gene"]
gtf.gene.plasmid <- gtf0.plasmid[gtf0.plasmid$type == "gene"]

head(gtf.gene.strains)
tail(gtf.gene.plasmid)

rm(gtf0.plasmid, gtf0.strains)


# Import VCF #### 

# name.plasmid <- bgzip("data/plasmid.filtered.vcf")
# name.strains <- bgzip("data/strains.filtered.vcf")

name.plasmid <- "data/plasmid.filtered.vcf.bgz"
name.strains <- "data/strains.filtered.vcf.bgz"

# indexTabix(name.plasmid, format = "vcf")
# indexTabix(name.strains, format = "vcf")

vcf.plasmid <- VcfFile(name.plasmid)
vcf.strains <- VcfFile(name.strains)

hdr.plasmid <- scanVcfHeader(vcf.plasmid)
hdr.strains <- scanVcfHeader(vcf.strains)

hdr.plasmid
hdr.strains

all.sam.plasmid <- samples(hdr.plasmid)
all.sam.strains <- samples(hdr.strains)

meta(hdr.plasmid)

for(i in names(meta(hdr.plasmid))){
  meta(hdr.plasmid)[[i]] %>% show()
}

info(hdr.plasmid)
info(hdr.strains)

info.desc <- info(hdr.plasmid)$Description
names(info.desc) <- rownames(info(hdr.plasmid))
View(info.desc)

readvcf.plasmid <- readVcf(vcf.plasmid)
readvcf.strains <- readVcf(vcf.strains)

save.image("01_dataimport.Rdata")

sessionInfo()
# R version 4.0.3 (2020-10-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
# [1] scater_1.18.3               SingleCellExperiment_1.12.0 ggplot2_3.3.2               snpStats_1.40.0             Matrix_1.2-18              
# [6] survival_3.2-7              VariantAnnotation_1.36.0    magrittr_2.0.1              SummarizedExperiment_1.20.0 MatrixGenerics_1.2.0       
# [11] matrixStats_0.57.0          Rsamtools_2.6.0             Biostrings_2.58.0           XVector_0.30.0              GenomicFeatures_1.42.1     
# [16] AnnotationDbi_1.52.0        Biobase_2.50.0              GenomicRanges_1.42.0        GenomeInfoDb_1.26.1         IRanges_2.24.0             
# [21] S4Vectors_0.28.0            BiocGenerics_0.36.0        
# 
# loaded via a namespace (and not attached):
#   [1] bitops_1.0-6              bit64_4.0.5               progress_1.2.2            httr_1.4.2                tools_4.0.3              
# [6] R6_2.5.0                  irlba_2.3.3               vipor_0.4.5               DBI_1.1.0                 colorspace_2.0-0         
# [11] withr_2.3.0               tidyselect_1.1.0          gridExtra_2.3             prettyunits_1.1.1         bit_4.0.4                
# [16] curl_4.3                  compiler_4.0.3            BiocNeighbors_1.8.1       xml2_1.3.2                DelayedArray_0.16.0      
# [21] rtracklayer_1.49.5        scales_1.1.1              askpass_1.1               rappdirs_0.3.1            stringr_1.4.0            
# [26] digest_0.6.27             pkgconfig_2.0.3           sparseMatrixStats_1.2.0   dbplyr_2.0.0              BSgenome_1.58.0          
# [31] rlang_0.4.9               rstudioapi_0.13           RSQLite_2.2.1             DelayedMatrixStats_1.12.1 generics_0.1.0           
# [36] BiocParallel_1.24.1       dplyr_1.0.2               RCurl_1.98-1.2            BiocSingular_1.6.0        GenomeInfoDbData_1.2.4   
# [41] scuttle_1.0.3             Rcpp_1.0.5                ggbeeswarm_0.6.0          munsell_0.5.0             viridis_0.5.1            
# [46] lifecycle_0.2.0           stringi_1.5.3             zlibbioc_1.36.0           BiocFileCache_1.14.0      grid_4.0.3               
# [51] blob_1.2.1                crayon_1.3.4              lattice_0.20-41           beachmat_2.6.2            splines_4.0.3            
# [56] hms_0.5.3                 pillar_1.4.7              biomaRt_2.46.0            XML_3.99-0.5              glue_1.4.2               
# [61] vctrs_0.3.5               gtable_0.3.0              openssl_1.4.3             purrr_0.3.4               assertthat_0.2.1         
# [66] xfun_0.19                 rsvd_1.0.3                viridisLite_0.3.0         tibble_3.0.4              GenomicAlignments_1.26.0 
# [71] tinytex_0.27              beeswarm_0.2.3            memoise_1.1.0             ellipsis_0.3.1           