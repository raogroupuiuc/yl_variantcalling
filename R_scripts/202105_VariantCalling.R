# Variant Calling analysis - Evolved Y. lipolytica PO1f mutants ####

# Redoing cleaner analysis with corrected reference genome May 03, 2021 

# Step 01 Set working directories and load libraries ####

setwd("D:/Box Sync/02 - RNA seq work/2101-YL-Evolution (Sangdo)/R analysis/")

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
library(ggpubr)

# Step 02a: Load data

# Load this datafilefor manipulation or go to step2 to run everything. 

load("202105_VariantCalling.Rdata")

# Step 02: Import and clean up genome information ####

rm(list = ls())
ref.strains <- "data/GCA_009372015.1_YarliW29_genomic.fna"
genome.strains <- FaFile(ref.strains)

ref.plasmid <- "data/PO1f_plasmid.fna"
genome.plasmid <- FaFile(ref.plasmid)

gtf0.strains <- rtracklayer::import("data/GCA_009372015.1_YarliW29_genomic.gff")
gtf0.plasmid <- rtracklayer::import("data/PO1f_plasmid.gff3")

gtf.gene.strains <- gtf0.strains[gtf0.strains$type == "gene"]
gtf.gene.plasmid <- gtf0.plasmid[gtf0.plasmid$type == "gene"]

rm(gtf0.plasmid, gtf0.strains)

# Step 03: Import VCF files and filter ####

name.plasmid <- "data/plasmid.filtered.vcf.bgz"
name.strains <- "data/strains.filtered.vcf.bgz"
# name.plasmid <- bgzip("data/plasmid.filtered.vcf")
# name.strains <- bgzip("data/strains.filtered.vcf")

# indexTabix(name.plasmid, format = "vcf")
# indexTabix(name.strains, format = "vcf")

vcf.plasmid <- VcfFile(name.plasmid)
vcf.strains <- VcfFile(name.strains)

hdr.plasmid <- scanVcfHeader(vcf.plasmid)

info.desc <- info(hdr.plasmid)$Description
names(info.desc) <- rownames(info(hdr.plasmid))
View(info.desc)

for (i in names(meta(hdr.plasmid))){
  meta(hdr.plasmid)[[i]] %>% show()
  }


readvcf.plasmid <- readVcf(vcf.plasmid)
readvcf.strains <- readVcf(vcf.strains)

# Filtering data for only PASSED values 

head(fixed(readvcf.plasmid))
# DataFrame with 6 rows and 4 columns
# REF                ALT      QUAL      FILTER
# <DNAStringSet> <DNAStringSetList> <numeric> <character>
#   1        ACCCCCCCCCCC                  A   2959.36   SORfilter
# 2      TAAAAAAAAAAAAA                  T   3912.36        PASS
# 3                   G                 GT    946.34        PASS
# 4                   T                 TC   2291.36   SORfilter
# 5 GAAAAAAAAAAAAAAAAAA                  G  12588.36   SORfilter
# 6            AGGGGGGG                  A   1938.18   SORfilter

head(fixed(readvcf.plasmid)$FILTER == "PASS")
# FALSE  TRUE  TRUE FALSE FALSE FALSE
(fixed(readvcf.plasmid)$FILTER == "PASS") %>% sum()
# 386
fixed(readvcf.plasmid)$FILTER %>% length()
# 815
(fixed(readvcf.strains)$FILTER == "PASS") %>% sum()
# 439
fixed(readvcf.strains)$FILTER %>% length()
# 879


# create a filter function for PASS 

pass01 <- function(vcf){
  fixed(vcf)$FILTER == "PASS"
}

filterVcf(vcf.plasmid, genome = "YL_PO1f_pINT03XYL123",
          destination = "data/plasmid_R_fil.vcf",
          filters = FilterRules(list(pass01)),
          index = TRUE)

filterVcf(vcf.strains, genome = "YL_PO1f",
          destination = "data/strains_R_fil.vcf",
          filters = FilterRules(list(pass01)),
          index = TRUE)

# Step 3b: Import filtered vcf files #### 


name.plasmid <- "data/plasmid.filtered.vcf.bgz"
name.strains <- "data/strains.filtered.vcf.bgz"
# name.plasmid <- bgzip("data/plasmid.filtered.vcf")
# name.strains <- bgzip("data/strains.filtered.vcf")

# indexTabix(name.plasmid, format = "vcf")
# indexTabix(name.strains, format = "vcf")

filt.vcf.plasmid <- VcfFile("data/plasmid_R_fil.vcf.bgz")
filt.vcf.strains <- VcfFile("data/strains_R_fil.vcf.bgz")

filt.readvcf.plasmid <- readVcf(filt.vcf.plasmid)
filt.readvcf.strains <- readVcf(filt.vcf.strains)


# Step 3c: Checking the filter cutoffs  #### 

# creating some plots  for checking filtering 

'Filtered using: 
https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants'

'Code: gatk VariantFiltration \
        -V plasmid.vcf \
        -O plasmid.filtered.vcf \
        -R ../../genome/PO1f_plasmid.fna \
        --cluster-size 3 \
        --cluster-window-size 10 \
        --filter-expression "QD < 2.0" --filter-name "QDFilter" \
        --filter-expression "MQ < 40.0" --filter-name "MQFilter" \
        --filter-expression "FS > 60.0" --filter-name "FSFilter" \
        --filter-expression "SOR > 3.0" --filter-name "SORfilter" \
        --filter-expression "HaplotypeScore > 13.0" --filter-name "HaplotypeScoreFilter" \
        --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSumFilter" \
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSumFilter" '

df.plt <- as.data.frame(info(readvcf.plasmid))
df.filt.plt <-as.data.frame(info(filt.readvcf.plasmid))
head(df.plt)
head(df.filt.plt)

# Qual By Depth 
' This is the variant confidence (from the QUAL field) divided by the unfiltered 
depth of non-hom-ref samples. This annotation is intended to normalize the 
variant quality in order to avoid inflation caused when there is deep coverage.
For filtering purposes it is better to use QD than either QUAL or DP directly.'


QD.plt <- ggplot(df.plt, aes(x = QD)) +
  geom_density(size = 1) + xlim(-10,50) + ggtitle("Prefilter") + 
  geom_vline(xintercept = 2, linetype="dashed", color = "#4682B4", size=0.75)

QD.plt.filt <- ggplot(df.filt.plt, aes(x = QD)) +
  geom_density(size = 1) + xlim(-10,50) + ggtitle("Filtered, 2") +
  geom_vline(xintercept = 2, linetype="dashed", color = "#4682B4", size=0.75)

a <- ggarrange(QD.plt, QD.plt.filt, ncol = 2, nrow = 1)

# RMSMappingQuality (MQ)
' This is the root mean square mapping quality over all the reads at the site. 
Instead of the average mapping quality of the site, this annotation gives the 
square root of the average of the squares of the mapping qualities at the site. 
It is meant to include the standard deviation of the mapping qualities. 
Including the standard deviation allows us to include the variation in the 
dataset. A low standard deviation means the values are all close to the mean, 
whereas a high standard deviation means the values are all far from the mean. 
When the mapping qualities are good at a site, the MQ will be around 60.'


MQ.plt <- ggplot(df.plt, aes(x = MQ)) +
  geom_density(size = 1) + ggtitle("Prefilter") + xlim(30,90) + 
  geom_vline(xintercept = 40, linetype="dashed", color = "#4682B4", size=0.75)

MQ.plt.filt <- ggplot(df.filt.plt, aes(x = MQ)) +
  geom_density(size = 1) + ggtitle("Filtered, 40") + xlim(30,90) +
  geom_vline(xintercept = 40, linetype="dashed", color = "#4682B4", size=0.75)

b <- ggarrange(MQ.plt, MQ.plt.filt, ncol = 2, nrow = 1)

# FisherStrand (FS)
' This is the Phred-scaled probability that there is strand bias at the site. 
Strand Bias tells us whether the alternate allele was seen more or less often 
on the forward or reverse strand than the reference allele. When there little 
to no strand bias at the site, the FS value will be close to 0.'


FS.plt <- ggplot(df.plt, aes(x = FS)) + scale_x_continuous(trans='log10') +
  geom_density(size = 1) + ggtitle("Prefilter") + #xlim(0.1,1000) + 
  geom_vline(xintercept = 60, linetype="dashed", color = "#4682B4", size=0.75)

FS.plt.filt <- ggplot(df.filt.plt, aes(x = FS)) + scale_x_continuous(trans='log10') +
  geom_density(size = 1) + ggtitle("Filtered, 60") + #xlim(0.1,1000) +
  geom_vline(xintercept = 60, linetype="dashed", color = "#4682B4", size=0.75)

c <- ggarrange(FS.plt, FS.plt.filt, ncol = 2, nrow = 1)


# StrandOddsRatio (SOR)
' This is another way to estimate strand bias using a test similar to the 
symmetric odds ratio test. SOR was created because FS tends to penalize variants 
that occur at the ends of exons. Reads at the ends of exons tend to only be 
covered by reads in one direction and FS gives those variants a bad score. SOR 
will take into account the ratios of reads that cover both alleles.'


SOR.plt <- ggplot(df.plt, aes(x = SOR)) +
  geom_density(size = 1) + ggtitle("Prefilter") + xlim(0,12) + 
  geom_vline(xintercept = 3, linetype="dashed", color = "#4682B4", size=0.75)

SOR.plt.filt <- ggplot(df.filt.plt, aes(x = SOR)) +
  geom_density(size = 1) + ggtitle("Filtered, 3") + xlim(0,12) + 
  geom_vline(xintercept = 3, linetype="dashed", color = "#4682B4", size=0.75)

d<- ggarrange(SOR.plt, SOR.plt.filt, ncol = 2, nrow = 1)

# MappingQualityRankSumTest (MQRankSum)
' This is the u-based z-approximation from the Rank Sum Test for mapping 
qualities. It compares the mapping qualities of the reads supporting the 
reference allele and the alternate allele. A positive value means the mapping 
qualities of the reads supporting the alternate allele are higher than those 
supporting the reference allele; a negative value indicates the mapping 
qualities of the reference allele are higher than those supporting the
alternate allele. A value close to zero is best and indicates little difference 
between the mapping qualities.'


MQRankSum.plt <- ggplot(df.plt, aes(x = MQRankSum)) +
  geom_density(size = 1) + ggtitle("Prefilter") + xlim(-50,50) + 
  geom_vline(xintercept = -12.5, linetype="dashed", color = "#4682B4", size=0.75)

MQRankSum.plt.filt <- ggplot(df.filt.plt, aes(x = MQRankSum)) +
  geom_density(size = 1) + ggtitle("Filtered, -12.5") + xlim(-50,50) + 
  geom_vline(xintercept = -12.5, linetype="dashed", color = "#4682B4", size=0.75)

e <- ggarrange(MQRankSum.plt, MQRankSum.plt.filt, ncol = 2, nrow = 1)

# ReadPosRankSumTest (ReadPosRankSum)
' The last annotation we will look at is ReadPosRankSum. This is the u-based 
z-approximation from the Rank Sum Test for site position within reads. It 
compares whether the positions of the reference and alternate alleles are 
different within the reads. Seeing an allele only near the ends of reads is 
indicative of error, because that is where sequencers tend to make the most 
errors. A negative value indicates that the alternate allele is found at the 
ends of reads more often than the reference allele; a positive value indicates 
that the reference allele is found at the ends of reads more often than the 
alternate allele. A value close to zero is best because it indicates there is 
little difference between the positions of the reference and alternate alleles 
in the reads. '


# x11(w=6, h=3)

ReadPosRankSum.plt <- ggplot(df.plt, aes(x = ReadPosRankSum)) +
  geom_density(size = 1) + ggtitle("Prefilter") + xlim(-20,30) + 
  geom_vline(xintercept = -8, linetype="dashed", color = "#4682B4", size=0.75)

ReadPosRankSum.plt.filt <- ggplot(df.filt.plt, aes(x = ReadPosRankSum)) +
  geom_density(size = 1) + ggtitle("Filtered, -8") + xlim(-20,30) + 
  geom_vline(xintercept = -8, linetype="dashed", color = "#4682B4", size=0.75)

f <- ggarrange(ReadPosRankSum.plt, ReadPosRankSum.plt.filt, ncol = 2, nrow = 1)

# x11(w=10, h=8)
jpeg("results/Filtercutoffs.jpeg", width = 10, height = 8, units = "in", res = 300)
ggarrange(a,b,c,d,e,f, ncol = 2, nrow = 3, 
          labels = c("A", "B", "C", "D", "E", "F"))
dev.off()

# Step 04: Identifying functional consequences of SNPs in CN evolved strains ####

txdb.strains <- makeTxDbFromGFF("data/GCA_009372015.1_YarliW29_genomic.gff",
                                format = "gff", 
                                organism = "Yarrowia lipolytica",
                                dataSource = "NCBI") # works

txdb.plasmid <- makeTxDbFromGFF("data/PO1f_plasmid.gff3",
                                format = "gff", 
                                organism = "Yarrowia lipolytica",
                                dataSource = "NCBI") # works

pc.strains <- predictCoding(filt.readvcf.strains, txdb.strains, genome.strains)
table(pc.strains$CONSEQUENCE)
# frameshift      nonsense nonsynonymous    synonymous 
# 6             3            45            32 

df.strains <- data.frame(SNP = names(pc.strains),
                         Chromosome = as.vector(seqnames(pc.strains)),
                         Position = start(pc.strains),
                         Transcript = pc.strains$TXID,
                         Gene = pc.strains$GENEID,
                         Consequence = pc.strains$CONSEQUENCE,
                         # AA_position = unlist(pc.strains$PROTEINLOC),
                         Ref_Codon = pc.strains$REFCODON,
                         Var_Codon = pc.strains$VARCODON, 
                         Ref_AA = pc.strains$REFAA,
                         Var_AA = pc.strains$VARAA)

# write.csv(df.strains, file = "results/protein_coding_variants_strains.csv", row.names = FALSE)


# Step 05: Identifying functional consequences of SNPs in Xylose evolved strains #### 
pc.plasmid <- predictCoding(filt.readvcf.plasmid, txdb.plasmid, genome.plasmid)

table(pc.plasmid$CONSEQUENCE)
# frameshift      nonsense nonsynonymous    synonymous 
# 6             5            27            21
df.plasmid <- data.frame(SNP = names(pc.plasmid),
                         Chromosome = as.vector(seqnames(pc.plasmid)),
                         Position = start(pc.plasmid),
                         Transcript = pc.plasmid$TXID,
                         Gene = pc.plasmid$GENEID,
                         Consequence = pc.plasmid$CONSEQUENCE,
                         # AA_position = unlist(pc.plasmid$PROTEINLOC),
                         Ref_Codon = pc.plasmid$REFCODON,
                         Var_Codon = pc.plasmid$VARCODON, 
                         Ref_AA = pc.plasmid$REFAA,
                         Var_AA = pc.plasmid$VARAA)
# write.csv(df.plasmid, file = "results/protein_coding_variants_plasmid.csv", row.names = FALSE)


# Step 06: List of mutations in each sample ####

gt.strains <- geno(filt.readvcf.strains)$GT
gt.plasmid <- geno(filt.readvcf.plasmid)$GT

table(gt.strains)
# .    0    1    2    3    4 
# 83  316 1722  69    4    1
table(gt.plasmid)
# .   0   1   2 
# 33  149 954  22 

colnames(gt.strains)
colnames(gt.plasmid)

mean(apply(gt.strains, 1, function(x) any(x %in% c("2", "3", "4")))) 
# 7.1% of loci have multiple alleles
mean(apply(gt.plasmid, 1, function(x) any(x %in% c("2")))) 
# 4.1% of loci have multiple alleles

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

mat.strains <- haploidSnpMat(filt.readvcf.strains)
mat.plasmid <- haploidSnpMat(filt.readvcf.plasmid)

summary(mat.plasmid)
summary(mat.strains)

# write.table(mat.plasmid, "results/snps-plasmid.txt", sep = "\t")
# write.table(mat.strains, "results/snps-strains.txt", sep = "\t")



# Step 07: Put everything together and output  ####

mat.plasmid <- as.data.frame(mat.plasmid)
mat.plasmid$SNP <- rownames(mat.plasmid)
final.plasmid <- merge(mat.plasmid, df.plasmid, by = "SNP", all = TRUE, sort = FALSE)

mat.strains <- as.data.frame(mat.strains)
mat.strains$SNP <- rownames(mat.strains)
final.strains <- merge(mat.strains, df.strains, by = "SNP", all = TRUE, sort = FALSE)

write.csv(final.plasmid, file = "results/mutations_plasmid.csv", row.names = FALSE)
write.csv(final.strains, file = "results/mutations_strains.csv", row.names = FALSE)

save.image("202105_VariantCalling.Rdata")

sessionInfo()
# R version 4.0.5 (2021-03-31)
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
#   [1] ggpubr_0.4.0                scater_1.18.6               SingleCellExperiment_1.12.0 ggplot2_3.3.3               snpStats_1.40.0            
# [6] Matrix_1.3-2                survival_3.2-10             VariantAnnotation_1.36.0    magrittr_2.0.1              SummarizedExperiment_1.20.0
# [11] MatrixGenerics_1.2.1        matrixStats_0.58.0          Rsamtools_2.6.0             Biostrings_2.58.0           XVector_0.30.0             
# [16] GenomicFeatures_1.42.2      AnnotationDbi_1.52.0        Biobase_2.50.0              GenomicRanges_1.42.0        GenomeInfoDb_1.26.4        
# [21] IRanges_2.24.1              S4Vectors_0.28.1            BiocGenerics_0.36.0        
# 
# loaded via a namespace (and not attached):
#   [1] ggbeeswarm_0.6.0          colorspace_2.0-0          ggsignif_0.6.1            ellipsis_0.3.1            rio_0.5.26               
# [6] scuttle_1.0.4             BiocNeighbors_1.8.2       rstudioapi_0.13           farver_2.1.0              bit64_4.0.5              
# [11] fansi_0.4.2               xml2_1.3.2                splines_4.0.5             sparseMatrixStats_1.2.1   cachem_1.0.4             
# [16] broom_0.7.6               dbplyr_2.1.0              compiler_4.0.5            httr_1.4.2                backports_1.2.1          
# [21] assertthat_0.2.1          fastmap_1.1.0             BiocSingular_1.6.0        prettyunits_1.1.1         tools_4.0.5              
# [26] rsvd_1.0.3                gtable_0.3.0              glue_1.4.2                GenomeInfoDbData_1.2.4    dplyr_1.0.5              
# [31] rappdirs_0.3.3            tinytex_0.31              Rcpp_1.0.6                carData_3.0-4             cellranger_1.1.0         
# [36] vctrs_0.3.6               rtracklayer_1.49.5        DelayedMatrixStats_1.12.3 xfun_0.22                 stringr_1.4.0            
# [41] openxlsx_4.2.3            beachmat_2.6.4            lifecycle_1.0.0           irlba_2.3.3               rstatix_0.7.0            
# [46] XML_3.99-0.6              zlibbioc_1.36.0           scales_1.1.1              BSgenome_1.58.0           hms_1.0.0                
# [51] curl_4.3                  memoise_2.0.0             gridExtra_2.3             biomaRt_2.46.3            stringi_1.5.3            
# [56] RSQLite_2.2.5             zip_2.1.1                 BiocParallel_1.24.1       rlang_0.4.10              pkgconfig_2.0.3          
# [61] bitops_1.0-6              lattice_0.20-41           purrr_0.3.4               GenomicAlignments_1.26.0  labeling_0.4.2           
# [66] cowplot_1.1.1             bit_4.0.4                 tidyselect_1.1.0          R6_2.5.0                  generics_0.1.0           
# [71] DelayedArray_0.16.3       DBI_1.1.1                 pillar_1.5.1              haven_2.4.1               foreign_0.8-81           
# [76] withr_2.4.1               abind_1.4-5               RCurl_1.98-1.3            tibble_3.1.0              crayon_1.4.1             
# [81] car_3.0-10                utf8_1.2.1                BiocFileCache_1.14.0      viridis_0.5.1             progress_1.2.2           
# [86] grid_4.0.5                readxl_1.3.1              data.table_1.14.0         blob_1.2.1                forcats_0.5.1            
# [91] digest_0.6.27             tidyr_1.1.3               openssl_1.4.3             munsell_0.5.0             beeswarm_0.3.1           
# [96] viridisLite_0.3.0         vipor_0.4.5               askpass_1.1 