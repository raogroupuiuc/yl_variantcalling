' Analysis of actual data 
Started March 4, 2021 from workshop practised code 
used for Xylose mutants 
Step 2- analyse the imported vcf files '

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


write.table(geno(readvcf.plasmid)$DP, "results/readdepth-plasmid.txt", sep = "\t")
write.table(geno(readvcf.strains)$DP, "results/readdepth-strains.txt", sep = "\t")

write.table(geno(readvcf.plasmid)$GT, "results/mutations-plasmid.txt", sep = "\t")
write.table(geno(readvcf.strains)$GT, "results/mutations-strains.txt", sep = "\t")


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

write.table(hits.plasmid4, "results/snpsingene-annotated-plasmid.txt", sep = "\t")
write.table(hits.strains4, "results/snpsingene-annotated-strains.txt", sep = "\t")


save.image("02-analysis-part1.Rdata")

# Identifying functional consequences of SNPs #### 

snps2search.plasmid <- subsetByOverlaps(sig.hits.plasmid, gtf.gene.plasmid)
snps2search.strains <- subsetByOverlaps(sig.hits.strains, gtf.gene.strains)

# Errors from here #### 
'Code after this does not work '

# genome.plasmid <- FaFile("data/pINT03-XYL123.fa")
txdb.strains <- makeTxDbFromGFF("data/GCA_009372015.1_YarliW29_genomic.gff",
                           format = "gff", organism = "Yarrowia lipolytica",
                           dataSource = "NCBI")

txdb.plasmid <- makeTxDbFromGFF("data/pINT03-XYL123.gff",
                                format = "gff", organism = "Yarrowia lipolytica",
                                dataSource = "NCBI")


genes(txdb.strains)

exonsBy(txdb.strains, by = "tx")

pc <- predictCoding(snps2search.strains, txdb.strains, genome.strains)
pc

table(pc$CONSEQUENCE)
my_txnames <- select(my_txdb, pc$TXID, c("TXID", "TXNAME"), "TXID")
head(my_txnames)
identical(my_txnames$TXID, as.integer(pc$TXID))
pc$TXNAME <- my_txnames$TXNAME

df <- data.frame(SNP = names(pc),
                 Chromosome = seqnames(pc),
                 Position = start(pc),
                 Transcript = pc$TXNAME,
                 Consequence = pc$CONSEQUENCE,
                 AA_position = unlist(pc$PROTEINLOC),
                 Ref_AA = pc$REFAA,
                 Var_AA = pc$VARAA)
write.csv(df, file = "protein_coding_variants.csv", row.names = FALSE)

# Converting from VCF to snpMatrix #### 

snpmat.plasmid <- genotypeToSnpMatrix(readvcf.plasmid)
snpmat.plasmid
# $genotypes
# A SnpMatrix with  3 rows and  171 columns
# Row names:  GEV_ACGCACCT-CCTTCACC ... XEV_GTATGTTC-TTCCTGTT 
# Col names:  ML755977.1:177460_ACCCCCCCCCCC/A ... ML756003.1:16734_G/GT 
# 
# $map
# DataFrame with 171 rows and 4 columns
# snp.names                allele.1           allele.2    ignore
# <character>          <DNAStringSet> <DNAStringSetList> <logical>
#   1   ML755977.1:177460_AC..            ACCCCCCCCCCC                  A      TRUE
# 2   ML755977.1:237488_TA..          TAAAAAAAAAAAAA                  T      TRUE
# 3   ML755977.1:351610_G/GT                       G                 GT      TRUE
# 4   ML755977.1:601795_T/TC                       T                 TC      TRUE
# 5   ML755978.1:83_GAAAAA..     GAAAAAAAAAAAAAAAAAA                  G      TRUE
# ...                    ...                     ...                ...       ...
# 167 ML756001.1:91052_CTT.. CTTTTTTTTT...TTTTTTTTTT                  C      TRUE
# 168 ML756001.1:93548_GAA.. GAAAAAAAAA...AAACAAAAAA                  G      TRUE
# 169 ML756001.1:196467_G/..                       G                GTT      TRUE
# 170 ML756001.1:235553_GG.. GGTTTTTTTT...GTTTCTTGTA                  G      TRUE
# 171  ML756003.1:16734_G/GT                       G             GT,GTT      TRUE

snpmat.strains <- genotypeToSnpMatrix(readvcf.strains)
snpmat.strains
# $genotypes
# A SnpMatrix with  5 rows and  143 columns
# Row names:  16EV_TCGATATC-ACTCGTGT ... PO1f_TATCGCAC-ACACTAAG 
# Col names:  ML755977.1:177460_ACCCCCCCCCCC/A ... ML755994.1:151700_GAAAA/G 
# 
# $map
# DataFrame with 143 rows and 4 columns
# snp.names        allele.1           allele.2    ignore
# <character>  <DNAStringSet> <DNAStringSetList> <logical>
#   1   ML755977.1:177460_AC..    ACCCCCCCCCCC                  A      TRUE
# 2    ML755977.1:212109_C/T               C                  T     FALSE
# 3   ML755977.1:237488_TA.. TAAAAAAAAAAAAAA               T,TA      TRUE
# 4   ML755977.1:351610_G/GT               G             GT,GTT      TRUE
# 5   ML755977.1:601795_T/TC               T                 TC      TRUE
# ...                    ...             ...                ...       ...
# 139  ML755994.1:34211_T/TC               T                 TC      TRUE
# 140  ML755994.1:42496_A/AC               A                 AC      TRUE
# 141  ML755994.1:139618_A/C               A                  C     FALSE
# 142  ML755994.1:139619_A/G               A                  G     FALSE
# 143 ML755994.1:151700_GA..           GAAAA              G,GAA      TRUE


as.numeric(snpmat.plasmid$genotypes)
as.numeric(snpmat.strains$genotypes)
# non-single nucleotide variations are set to NA
# Warning messages:
#   1: In .local(x, ...) : variants with >1 ALT allele are set to NA
# 2: In .local(x, ...) : non-diploid variants are set to NA

# cannot do statistics and filtering because of ignore filter in genotypetoSnpMatrix 

# Basic statistics and filtering on SNPs #### 

mat <- mysnpmat$genotypes[,!mysnpmat$map$ignore]
mat

summary(mat)

sample_stats <- row.summary(mat)

ggplot(sample_stats, aes(x = Call.rate, y = Heterozygosity)) +
  geom_point()

# filtering by heterozygosity 

highhet <- isOutlier(sample_stats$Heterozygosity, type = "higher")

ggplot(sample_stats, aes(x = Call.rate, y = Heterozygosity,
                         color = highhet)) +
  geom_point()

mat2 <- mat[!highhet,]
mat2

marker_stats <- col.summary(mat2)

ggplot(marker_stats, aes(x = MAF)) +
  geom_histogram()

marker_stats <- dplyr::mutate(marker_stats, Ho.He = P.AB / (2 * MAF * (1 - MAF)))

ggplot(marker_stats, aes(x = Ho.He)) +
  geom_histogram()

highhet2 <- isOutlier(marker_stats$Ho.He, type = "higher")

highhet2[is.na(highhet2)] <- TRUE

ggplot(marker_stats, aes(x = RAF, y = P.BB, color = highhet2)) + geom_point()

mat3 <- mat2[, !highhet2]
mat3

# Linkage Disequilibrium #### 

mydepth <- 100 # how many adjacent markers to look at
myLD <- ld(mat3, depth = mydepth, stats = "R.squared", symmetric = FALSE)

image(myLD[1:500, 1:500], lwd = 0)

pos <- start(rowRanges(mydata)[colnames(mat3)])
nSNP <- length(pos)
diags <- vector("list", mydepth)
for (i in 1:mydepth) diags[[i]] <- pos[(i+1):nSNP] - pos[1:(nSNP-i)]
physical_distance <- bandSparse(nSNP, k=1:mydepth, diagonals=diags)

physical_distance_vals <- physical_distance@x
LD_vals <- myLD@x

random_subset <- sample(which(physical_distance_vals < 2e5), 5000)

ggplot(mapping = aes(x = physical_distance_vals[random_subset],
                     y = LD_vals[random_subset])) +
  labs(x = "Distance in bp", y = "R-squared") +
  geom_point(alpha = 0.1) +
  geom_smooth(formula = y ~ log(x)) +
  geom_vline(xintercept = c(2000, 5000), color = "red", lty = 2)


# Principal components analysis #### 

my_xxt <- xxt(mat3)
my_pca <- eigen(my_xxt, symmetric = TRUE)
percent_variation <- round(my_pca$values/sum(my_pca$values) * 100, 2)
plot(percent_variation)

plotPCs <- function(x, y, eigenvect = my_pca$vectors,
                    pct_var = percent_variation){
  ggplot(mapping = aes(x = eigenvect[,x], y = eigenvect[,y])) +
    geom_point() +
    labs(x = paste0("PC", x, " (", pct_var[x], "%)"),
         y = paste0("PC", y, " (", pct_var[y], "%)"))
}

plotPCs(1, 2)
plotPCs(3, 4)
plotPCs(5, 6)

pca_tab <- data.frame(Sample = rownames(mat3), my_pca$vectors[,1:6])
colnames(pca_tab)[-1] <- paste0("PC", 1:6)

pca_tab %>% dplyr::filter(PC1 < -0.05) %>%
  dplyr::select(Sample, PC1, PC2)

# ISsue CODE ends #### 

# Surplus code from workshop - unused #### 

# Specify coordinates in a genome #### 
myqtl <- GRanges(c("Chr2", "Chr2", "Chr8"),
                 IRanges(start = c(134620000, 48023000, 150341000),
                         end   = c(134752000, 48046000, 150372000)))
myqtl
# GRanges object with 3 ranges and 0 metadata columns:
#   seqnames              ranges strand
# <Rle>           <IRanges>  <Rle>
#   [1]     Chr2 134620000-134752000      *
#   [2]     Chr2   48023000-48046000      *
#   [3]     Chr8 150341000-150372000      *
#   seqinfo: 2 sequences from an unspecified genome; no seqlengths

names(myqtl) <- c("Yld1", "LA1", "LA2")
myqtl$Trait <- c("Yield", "Leaf angle", "Leaf angle")
myqtl
# GRanges object with 3 ranges and 1 metadata column:
#   seqnames              ranges strand |       Trait
# <Rle>           <IRanges>  <Rle> | <character>
#   Yld1     Chr2 134620000-134752000      * |       Yield
# LA1     Chr2   48023000-48046000      * |  Leaf angle
# LA2     Chr8 150341000-150372000      * |  Leaf angle
#   seqinfo: 2 sequences from an unspecified genome; no seqlengths

myqtl[1]
myqtl["LA2"]
myqtl[myqtl$Trait == "Leaf angle"]

# Finding overlaps #### 

qtl_genes <- subsetByOverlaps(gtf1, myqtl)
qtl_genes


# DNA Sequences #### 

qtl_genes_seq <- scanFa(mygenome, qtl_genes)
qtl_genes_seq

reverseComplement(qtl_genes_seq)

translate(qtl_genes_seq)

# Extracting transcript sequences ####

' If you do want to get the sequences of transcripts or CDS from a genome, \
see ?GenomicFeatures::extractTranscriptSeqs. You would need to import the GFF 
with makeTxDbFromGFF rather than rtracklayer::import. If you need to create a 
DNAStringSet from scratch, you can do it directly from a character vector. '

test_dna <- DNAStringSet(c("AGGG", "TCAGATTTAAC", "TC"))
test_dna

'How would you import the full sequence for chromosome 3?'

chr3length <- seqlengths(seqinfo(mygenome))["Chr3"]
chr3range <- GRanges("Chr3", IRanges(start = 1, end = chr3length))
chr3seq <- scanFa(mygenome, chr3range)
chr3seq

# Experimental results #### 

data(airway, package="airway")

airway
colData(airway)
assays(airway)$counts %>% head(10)
rowRanges(airway)

myregion <- GRanges("4", IRanges(start = 82660000, end = 119280000))
airway4 <- subsetByOverlaps(airway, myregion)
airway4
rowRanges(airway4)

