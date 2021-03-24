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
genome.strains
# class: FaFile 
# path: data/GCA_009372015.1_YarliW29_genomic.fna
# index: data/GCA_009372015.1_YarliW29_genomic.fna.fai
# gzindex: data/GCA_009372015.1_YarliW29_genomic.fna.gzi
# isOpen: FALSE 
# yieldSize: NA 
seqinfo(genome.strains)
# Seqinfo object with 168 sequences from an unspecified genome:
#   seqnames   seqlengths isCircular genome
# ML755977.1     780106       <NA>   <NA>
#   ML755978.1     515156       <NA>   <NA>
#   ML755979.1     507099       <NA>   <NA>
#   ML755980.1     467546       <NA>   <NA>
#   ML755981.1     441470       <NA>   <NA>
#   ...               ...        ...    ...
# ML756140.1       1305       <NA>   <NA>
#   ML756141.1       1284       <NA>   <NA>
#   ML756142.1       1215       <NA>   <NA>
#   ML756143.1       1089       <NA>   <NA>
#   ML756144.1       1085       <NA>   <NA>
seqnames(seqinfo(genome.strains))
# [1] "ML755977.1" "ML755978.1" "ML755979.1" "ML755980.1" "ML755981.1" "ML755982.1"
# [7] "ML755983.1" "ML755984.1" "ML755985.1" "ML755986.1" "ML755987.1" "ML755988.1"
# [13] "ML755989.1" "ML755990.1" "ML755991.1" "ML755992.1" "ML755993.1" "ML755994.1"
# [19] "ML755995.1" "ML755996.1" "ML755997.1" "ML755998.1" "ML755999.1" "ML756000.1"
# [25] "ML756001.1" "ML756002.1" "ML756003.1" "ML756004.1" "ML756005.1" "ML756006.1"
# [31] "ML756007.1" "ML756008.1" "ML756009.1" "ML756010.1" "ML756011.1" "ML756012.1"
# [37] "ML756013.1" "ML756014.1" "ML756015.1" "ML756016.1" "ML756017.1" "ML756018.1"
# [43] "ML756019.1" "ML756020.1" "ML756021.1" "ML756022.1" "ML756023.1" "ML756024.1"
# [49] "ML756025.1" "ML756026.1" "ML756027.1" "ML756028.1" "ML756029.1" "ML756030.1"
# [55] "ML756031.1" "ML756032.1" "ML756033.1" "ML756034.1" "ML756035.1" "ML756036.1"
# [61] "ML756037.1" "ML756038.1" "ML756039.1" "ML756040.1" "ML756041.1" "ML756042.1"
# [67] "ML756043.1" "ML756044.1" "ML756045.1" "ML756046.1" "ML756047.1" "ML756048.1"
# [73] "ML756049.1" "ML756050.1" "ML756051.1" "ML756052.1" "ML756053.1" "ML756054.1"
# [79] "ML756055.1" "ML756056.1" "ML756057.1" "ML756058.1" "ML756059.1" "ML756060.1"
# [85] "ML756061.1" "ML756062.1" "ML756063.1" "ML756064.1" "ML756065.1" "ML756066.1"
# [91] "ML756067.1" "ML756068.1" "ML756069.1" "ML756070.1" "ML756071.1" "ML756072.1"
# [97] "ML756073.1" "ML756074.1" "ML756075.1" "ML756076.1" "ML756077.1" "ML756078.1"
# [103] "ML756079.1" "ML756080.1" "ML756081.1" "ML756082.1" "ML756083.1" "ML756084.1"
# [109] "ML756085.1" "ML756086.1" "ML756087.1" "ML756088.1" "ML756089.1" "ML756090.1"
# [115] "ML756091.1" "ML756092.1" "ML756093.1" "ML756094.1" "ML756095.1" "ML756096.1"
# [121] "ML756097.1" "ML756098.1" "ML756099.1" "ML756100.1" "ML756101.1" "ML756102.1"
# [127] "ML756103.1" "ML756104.1" "ML756105.1" "ML756106.1" "ML756107.1" "ML756108.1"
# [133] "ML756109.1" "ML756110.1" "ML756111.1" "ML756112.1" "ML756113.1" "ML756114.1"
# [139] "ML756115.1" "ML756116.1" "ML756117.1" "ML756118.1" "ML756119.1" "ML756120.1"
# [145] "ML756121.1" "ML756122.1" "ML756123.1" "ML756124.1" "ML756125.1" "ML756126.1"
# [151] "ML756127.1" "ML756128.1" "ML756129.1" "ML756130.1" "ML756131.1" "ML756132.1"
# [157] "ML756133.1" "ML756134.1" "ML756135.1" "ML756136.1" "ML756137.1" "ML756138.1"
# [163] "ML756139.1" "ML756140.1" "ML756141.1" "ML756142.1" "ML756143.1" "ML756144.1"

ref.plasmid <- "data/PO1f_plasmid.fna"
# indexFa(ref.plasmid)
genome.plasmid <- FaFile(ref.plasmid)
genome.plasmid
# class: FaFile 
# path: data/PO1f_plasmid.fna
# index: data/PO1f_plasmid.fna.fai
# gzindex: data/PO1f_plasmid.fna.gzi
# isOpen: FALSE 
# yieldSize: NA 
seqinfo(genome.plasmid)
# Seqinfo object with 169 sequences from an unspecified genome:
#   seqnames         seqlengths isCircular genome
# ML755977.1           780106       <NA>   <NA>
#   ML755978.1           515156       <NA>   <NA>
#   ML755979.1           507099       <NA>   <NA>
#   ML755980.1           467546       <NA>   <NA>
#   ML755981.1           441470       <NA>   <NA>
#   ...                     ...        ...    ...
# ML756141.1             1284       <NA>   <NA>
#   ML756142.1             1215       <NA>   <NA>
#   ML756143.1             1089       <NA>   <NA>
#   ML756144.1             1085       <NA>   <NA>
#   pINT03-XYL123.gb      13224       <NA>   <NA>
seqnames(seqinfo(genome.plasmid))
# [1] "ML755977.1" "ML755978.1" "ML755979.1" "ML755980.1" "ML755981.1" "ML755982.1"
# [7] "ML755983.1" "ML755984.1" "ML755985.1" "ML755986.1" "ML755987.1" "ML755988.1"
# [13] "ML755989.1" "ML755990.1" "ML755991.1" "ML755992.1" "ML755993.1" "ML755994.1"
# [19] "ML755995.1" "ML755996.1" "ML755997.1" "ML755998.1" "ML755999.1" "ML756000.1"
# [25] "ML756001.1" "ML756002.1" "ML756003.1" "ML756004.1" "ML756005.1" "ML756006.1"
# [31] "ML756007.1" "ML756008.1" "ML756009.1" "ML756010.1" "ML756011.1" "ML756012.1"
# [37] "ML756013.1" "ML756014.1" "ML756015.1" "ML756016.1" "ML756017.1" "ML756018.1"
# [43] "ML756019.1" "ML756020.1" "ML756021.1" "ML756022.1" "ML756023.1" "ML756024.1"
# [49] "ML756025.1" "ML756026.1" "ML756027.1" "ML756028.1" "ML756029.1" "ML756030.1"
# [55] "ML756031.1" "ML756032.1" "ML756033.1" "ML756034.1" "ML756035.1" "ML756036.1"
# [61] "ML756037.1" "ML756038.1" "ML756039.1" "ML756040.1" "ML756041.1" "ML756042.1"
# [67] "ML756043.1" "ML756044.1" "ML756045.1" "ML756046.1" "ML756047.1" "ML756048.1"
# [73] "ML756049.1" "ML756050.1" "ML756051.1" "ML756052.1" "ML756053.1" "ML756054.1"
# [79] "ML756055.1" "ML756056.1" "ML756057.1" "ML756058.1" "ML756059.1" "ML756060.1"
# [85] "ML756061.1" "ML756062.1" "ML756063.1" "ML756064.1" "ML756065.1" "ML756066.1"
# [91] "ML756067.1" "ML756068.1" "ML756069.1" "ML756070.1" "ML756071.1" "ML756072.1"
# [97] "ML756073.1" "ML756074.1" "ML756075.1" "ML756076.1" "ML756077.1" "ML756078.1"
# [103] "ML756079.1" "ML756080.1" "ML756081.1" "ML756082.1" "ML756083.1" "ML756084.1"
# [109] "ML756085.1" "ML756086.1" "ML756087.1" "ML756088.1" "ML756089.1" "ML756090.1"
# [115] "ML756091.1" "ML756092.1" "ML756093.1" "ML756094.1" "ML756095.1" "ML756096.1"
# [121] "ML756097.1" "ML756098.1" "ML756099.1" "ML756100.1" "ML756101.1" "ML756102.1"
# [127] "ML756103.1" "ML756104.1" "ML756105.1" "ML756106.1" "ML756107.1" "ML756108.1"
# [133] "ML756109.1" "ML756110.1" "ML756111.1" "ML756112.1" "ML756113.1" "ML756114.1"
# [139] "ML756115.1" "ML756116.1" "ML756117.1" "ML756118.1" "ML756119.1" "ML756120.1"
# [145] "ML756121.1" "ML756122.1" "ML756123.1" "ML756124.1" "ML756125.1" "ML756126.1"
# [151] "ML756127.1" "ML756128.1" "ML756129.1" "ML756130.1" "ML756131.1" "ML756132.1"
# [157] "ML756133.1" "ML756134.1" "ML756135.1" "ML756136.1" "ML756137.1" "ML756138.1"
# [163] "ML756139.1" "ML756140.1" "ML756141.1" "ML756142.1" "ML756143.1" "ML756144.1"
# [169] "pINT03-XYL123.gb"



# Import gene annotations #### 

gtf0.strains <- rtracklayer::import("data/GCA_009372015.1_YarliW29_genomic.gff")
gtf0.plasmid <- rtracklayer::import("data/PO1f_plasmid.gff")

gtf.gene.strains <- gtf0.strains[gtf0.strains$type == "gene"]
gtf.gene.plasmid <- gtf0.plasmid[gtf0.plasmid$type == "gene"]

head(gtf.gene.strains)
# GRanges object with 6 ranges and 30 metadata columns:
#   seqnames     ranges strand |   source     type     score     phase                     ID          Dbxref              Name  chromosome
# <Rle>  <IRanges>  <Rle> | <factor> <factor> <numeric> <integer>            <character> <CharacterList>       <character> <character>
#   [1] ML755977.1     21-780      + |  Genbank     gene        NA      <NA> gene-BKA90DRAFT_143812                 BKA90DRAFT_143812        <NA>
#   [2] ML755977.1   900-2824      - |  Genbank     gene        NA      <NA> gene-BKA90DRAFT_132270                 BKA90DRAFT_132270        <NA>
#   [3] ML755977.1  3366-4966      + |  Genbank     gene        NA      <NA> gene-BKA90DRAFT_132271                 BKA90DRAFT_132271        <NA>
#   [4] ML755977.1  5060-6289      - |  Genbank     gene        NA      <NA> gene-BKA90DRAFT_132272                 BKA90DRAFT_132272        <NA>
#   [5] ML755977.1  7888-9231      + |  Genbank     gene        NA      <NA>     gene-BKA90DRAFT_19                     BKA90DRAFT_19        <NA>
#   [6] ML755977.1 9339-11270      - |  Genbank     gene        NA      <NA> gene-BKA90DRAFT_150057                 BKA90DRAFT_150057        <NA>
#   gbkey      genome    mol_type      strain       end_range   gene_biotype         locus_tag     partial     start_range
# <character> <character> <character> <character> <CharacterList>    <character>       <character> <character> <CharacterList>
#   [1]        Gene        <NA>        <NA>        <NA>           780,. protein_coding BKA90DRAFT_143812        true            .,21
# [2]        Gene        <NA>        <NA>        <NA>                 protein_coding BKA90DRAFT_132270        <NA>                
#   [3]        Gene        <NA>        <NA>        <NA>                 protein_coding BKA90DRAFT_132271        <NA>                
#   [4]        Gene        <NA>        <NA>        <NA>                 protein_coding BKA90DRAFT_132272        <NA>                
#   [5]        Gene        <NA>        <NA>        <NA>                 protein_coding     BKA90DRAFT_19        <NA>                
#   [6]        Gene        <NA>        <NA>        <NA>         11270,. protein_coding BKA90DRAFT_150057        true          .,9339
# Parent orig_protein_id orig_transcript_id     product  protein_id   Ontology_term go_function go_component  go_process
# <CharacterList>     <character>        <character> <character> <character> <CharacterList> <character>  <character> <character>
#   [1]                            <NA>               <NA>        <NA>        <NA>                        <NA>         <NA>        <NA>
#   [2]                            <NA>               <NA>        <NA>        <NA>                        <NA>         <NA>        <NA>
#   [3]                            <NA>               <NA>        <NA>        <NA>                        <NA>         <NA>        <NA>
#   [4]                            <NA>               <NA>        <NA>        <NA>                        <NA>         <NA>        <NA>
#   [5]                            <NA>               <NA>        <NA>        <NA>                        <NA>         <NA>        <NA>
#   [6]                            <NA>               <NA>        <NA>        <NA>                        <NA>         <NA>        <NA>
#   Note   anticodon   exception      pseudo
# <CharacterList> <character> <character> <character>
#   [1]                        <NA>        <NA>        <NA>
#   [2]                        <NA>        <NA>        <NA>
#   [3]                        <NA>        <NA>        <NA>
#   [4]                        <NA>        <NA>        <NA>
#   [5]                        <NA>        <NA>        <NA>
#   [6]                        <NA>        <NA>        <NA>
#   seqinfo: 168 sequences from an unspecified genome; no seqlengths

tail(gtf.gene.plasmid)
# GRanges object with 6 ranges and 30 metadata columns:
#   seqnames      ranges strand |   source     type     score     phase          ID          Dbxref        Name  chromosome
# <Rle>   <IRanges>  <Rle> | <factor> <factor> <numeric> <integer> <character> <CharacterList> <character> <character>
#   [1] pINT03-XYL123.gb   4079-5035      + |  Genbank     gene        NA         0        XYL1                        <NA>        <NA>
#   [2] pINT03-XYL123.gb   5336-6302      + |  Genbank     gene        NA         0       TDH1P                        <NA>        <NA>
#   [3] pINT03-XYL123.gb   6303-7394      + |  Genbank     gene        NA         0        XYL2                        <NA>        <NA>
#   [4] pINT03-XYL123.gb   7695-8526      + |  Genbank     gene        NA         0       FBA1p                        <NA>        <NA>
#   [5] pINT03-XYL123.gb  8527-10398      + |  Genbank     gene        NA         0        XYL3                        <NA>        <NA>
#   [6] pINT03-XYL123.gb 10399-10698      + |  Genbank     gene        NA         0       FBA1t                        <NA>        <NA>
#   gbkey      genome    mol_type      strain       end_range gene_biotype   locus_tag     partial     start_range          Parent
# <character> <character> <character> <character> <CharacterList>  <character> <character> <character> <CharacterList> <CharacterList>
#   [1]        <NA>        <NA>        <NA>        <NA>                         <NA>        <NA>        <NA>                                
#   [2]        <NA>        <NA>        <NA>        <NA>                         <NA>        <NA>        <NA>                                
#   [3]        <NA>        <NA>        <NA>        <NA>                         <NA>        <NA>        <NA>                                
#   [4]        <NA>        <NA>        <NA>        <NA>                         <NA>        <NA>        <NA>                                
#   [5]        <NA>        <NA>        <NA>        <NA>                         <NA>        <NA>        <NA>                                
#   [6]        <NA>        <NA>        <NA>        <NA>                         <NA>        <NA>        <NA>                                
#   orig_protein_id orig_transcript_id     product  protein_id   Ontology_term go_function go_component  go_process            Note
# <character>        <character> <character> <character> <CharacterList> <character>  <character> <character> <CharacterList>
#   [1]            <NA>               <NA>        <NA>        <NA>                        <NA>         <NA>        <NA>                
#   [2]            <NA>               <NA>        <NA>        <NA>                        <NA>         <NA>        <NA>                
#   [3]            <NA>               <NA>        <NA>        <NA>                        <NA>         <NA>        <NA>                
#   [4]            <NA>               <NA>        <NA>        <NA>                        <NA>         <NA>        <NA>                
#   [5]            <NA>               <NA>        <NA>        <NA>                        <NA>         <NA>        <NA>                
#   [6]            <NA>               <NA>        <NA>        <NA>                        <NA>         <NA>        <NA>                
#   anticodon   exception      pseudo
# <character> <character> <character>
#   [1]        <NA>        <NA>        <NA>
#   [2]        <NA>        <NA>        <NA>
#   [3]        <NA>        <NA>        <NA>
#   [4]        <NA>        <NA>        <NA>
#   [5]        <NA>        <NA>        <NA>
#   [6]        <NA>        <NA>        <NA>
#   seqinfo: 169 sequences from an unspecified genome; no seqlengths

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
# class: VCFHeader 
# samples(3): GEV_ACGCACCT-CCTTCACC X123_CGCTATGT-GTGTCGGA XEV_GTATGTTC-TTCCTGTT
# meta(7): fileformat source ... GATKCommandLine contig
# fixed(2): FILTER ALT
# info(18): AC AF ... ReadPosRankSum SOR
# geno(8): GT AD ... RGQ SB

hdr.strains
# class: VCFHeader 
# samples(5): 16EV_TCGATATC-ACTCGTGT 1EV_TACTCATA-GCCACAGG 2EV_CGTCTGCG-ATTGTGAA 41EV_CTAGCGCT-GTCTACAC PO1f_TATCGCAC-ACACTAAG
# meta(7): fileformat source ... GATKCommandLine contig
# fixed(2): FILTER ALT
# info(18): AC AF ... ReadPosRankSum SOR
# geno(8): GT AD ... RGQ SB

all.sam.plasmid <- samples(hdr.plasmid)
all.sam.strains <- samples(hdr.strains)

meta(hdr.plasmid)
# DataFrameList of length 7
# names(7): fileformat source source.1 source.2 source.3 GATKCommandLine contig

for(i in names(meta(hdr.plasmid))){
  meta(hdr.plasmid)[[i]] %>% show()
}
# DataFrame with 1 row and 1 column
# Value
# <character>
#   fileformat     VCFv4.2
# DataFrame with 1 row and 1 column
# Value
# <character>
#   source CombineGVCFs
# DataFrame with 1 row and 1 column
# Value
# <character>
#   source.1 GenotypeGVCFs
# DataFrame with 1 row and 1 column
# Value
# <character>
#   source.2 HaplotypeCaller
# DataFrame with 1 row and 1 column
# Value
# <character>
#   source.3 VariantFiltration
# DataFrame with 4 rows and 3 columns
# CommandLine     Version                   Date
# <character> <character>            <character>
#   CombineGVCFs      "CombineGVCFs  --out..   "4.1.4.0" "February 25, 2021 4..
# GenotypeGVCFs     "GenotypeGVCFs  --ou..   "4.1.4.0" "February 25, 2021 4..
# HaplotypeCaller   "HaplotypeCaller  --..   "4.1.4.0" "February 25, 2021 1..
# VariantFiltration "VariantFiltration  ..   "4.1.4.0" "February 25, 2021 4..
# DataFrame with 168 rows and 1 column
# length
# <character>
#   ML755977.1      780106
# ML755978.1      515156
# ML755979.1      507099
# ML755980.1      467546
# ML755981.1      441470
# ...                ...
# ML756140.1        1305
# ML756141.1        1284
# ML756142.1        1215
# ML756143.1        1089
# ML756144.1        1085


info(hdr.plasmid)
info(hdr.strains)

info.desc <- info(hdr.plasmid)$Description
names(info.desc) <- rownames(info(hdr.plasmid))
View(info.desc)

# load("01_dataimport.Rdata")

readvcf.plasmid <- readVcf(vcf.plasmid)
# class: CollapsedVCF 
# dim: 171 3 
# rowRanges(vcf):
#   GRanges with 5 metadata columns: paramRangeID, REF, ALT, QUAL, FILTER
# info(vcf):
#   DataFrame with 18 columns: AC, AF, AN, BaseQRankSum, DP, DS, END, ExcessHet, FS, Inbreeding...
# info(header(vcf)):
#   Number Type    Description                                                   
# AC              A      Integer Allele count in genotypes, for each ALT allele, in the same...
# AF              A      Float   Allele Frequency, for each ALT allele, in the same order as...
# AN              1      Integer Total number of alleles in called genotypes                   
# BaseQRankSum    1      Float   Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qua...
# DP              1      Integer Approximate read depth; some reads may have been filtered     
# DS              0      Flag    Were any of the samples downsampled?                          
# END             1      Integer Stop position of the interval                                 
# ExcessHet       1      Float   Phred-scaled p-value for exact test of excess heterozygosity  
# FS              1      Float   Phred-scaled p-value using Fisher's exact test to detect st...
# InbreedingCoeff 1      Float   Inbreeding coefficient as estimated from the genotype likel...
# MLEAC           A      Integer Maximum likelihood expectation (MLE) for the allele counts ...
# MLEAF           A      Float   Maximum likelihood expectation (MLE) for the allele frequen...
# MQ              1      Float   RMS Mapping Quality                                           
# MQRankSum       1      Float   Z-score From Wilcoxon rank sum test of Alt vs. Ref read map...
# QD              1      Float   Variant Confidence/Quality by Depth                           
# RAW_MQandDP     2      Integer Raw data (sum of squared MQ and total depth) for improved R...
# ReadPosRankSum  1      Float   Z-score from Wilcoxon rank sum test of Alt vs. Ref read pos...
# SOR             1      Float   Symmetric Odds Ratio of 2x2 contingency table to detect str...
# geno(vcf):
#   List of length 8: GT, AD, DP, GQ, MIN_DP, PL, RGQ, SB
# geno(header(vcf)):
#           Number Type    Description                                                            
#    GT     1      String  Genotype                                                               
#    AD     R      Integer Allelic depths for the ref and alt alleles in the order listed         
#    DP     1      Integer Approximate read depth (reads with MQ=255 or with bad mates are filt...
#    GQ     1      Integer Genotype Quality                                                       
#    MIN_DP 1      Integer Minimum DP observed within the GVCF block                              
#    PL     G      Integer Normalized, Phred-scaled likelihoods for genotypes as defined in the...
#    RGQ    1      Integer Unconditional reference genotype confidence, encoded as a phred qual...
#    SB     4      Integer Per-sample component statistics which comprise the Fisher's Exact Te...


readvcf.strains <- readVcf(vcf.strains)
# class: CollapsedVCF 
# dim: 143 5 
# rowRanges(vcf):
#   GRanges with 5 metadata columns: paramRangeID, REF, ALT, QUAL, FILTER
# info(vcf):
#   DataFrame with 18 columns: AC, AF, AN, BaseQRankSum, DP, DS, END, ExcessHet, FS, Inbr...
# info(header(vcf)):
#   Number Type    Description                                             
# AC              A      Integer Allele count in genotypes, for each ALT allele, in th...
# AF              A      Float   Allele Frequency, for each ALT allele, in the same or...
# AN              1      Integer Total number of alleles in called genotypes             
# BaseQRankSum    1      Float   Z-score from Wilcoxon rank sum test of Alt Vs. Ref ba...
# DP              1      Integer Approximate read depth; some reads may have been filt...
# DS              0      Flag    Were any of the samples downsampled?                    
#   END             1      Integer Stop position of the interval                           
# ExcessHet       1      Float   Phred-scaled p-value for exact test of excess heteroz...
# FS              1      Float   Phred-scaled p-value using Fisher's exact test to det...
#    InbreedingCoeff 1      Float   Inbreeding coefficient as estimated from the genotype...
#    MLEAC           A      Integer Maximum likelihood expectation (MLE) for the allele c...
#    MLEAF           A      Float   Maximum likelihood expectation (MLE) for the allele f...
#    MQ              1      Float   RMS Mapping Quality                                     
#    MQRankSum       1      Float   Z-score From Wilcoxon rank sum test of Alt vs. Ref re...
#    QD              1      Float   Variant Confidence/Quality by Depth                     
#    RAW_MQandDP     2      Integer Raw data (sum of squared MQ and total depth) for impr...
#    ReadPosRankSum  1      Float   Z-score from Wilcoxon rank sum test of Alt vs. Ref re...
#    SOR             1      Float   Symmetric Odds Ratio of 2x2 contingency table to dete...
# geno(vcf):
#   List of length 8: GT, AD, DP, GQ, MIN_DP, PL, RGQ, SB
# geno(header(vcf)):
#           Number Type    Description                                                      
#    GT     1      String  Genotype                                                         
#    AD     R      Integer Allelic depths for the ref and alt alleles in the order listed   
#    DP     1      Integer Approximate read depth (reads with MQ=255 or with bad mates ar...
#    GQ     1      Integer Genotype Quality                                                 
#    MIN_DP 1      Integer Minimum DP observed within the GVCF block                        
#    PL     G      Integer Normalized, Phred-scaled likelihoods for genotypes as defined ...
#    RGQ    1      Integer Unconditional reference genotype confidence, encoded as a phre...
#    SB     4      Integer Per-sample component statistics which comprise the Fisher's Ex...

head(rowRanges(readvcf.plasmid))
# GRanges object with 6 ranges and 5 metadata columns:
#   seqnames        ranges strand | paramRangeID
# <Rle>     <IRanges>  <Rle> |     <factor>
#   ML755977.1:177460_ACCCCCCCCCCC/A ML755977.1 177460-177471      * |           NA
# ML755977.1:237488_TAAAAAAAAAAAAA/T ML755977.1 237488-237501      * |           NA
# ML755977.1:351610_G/GT ML755977.1        351610      * |           NA
# ML755977.1:601795_T/TC ML755977.1        601795      * |           NA
# ML755978.1:83_GAAAAAAAAAAAAAAAAAA/G ML755978.1        83-101      * |           NA
# ML755978.1:66627_AGGGGGGG/A ML755978.1   66627-66634      * |           NA
# REF                ALT      QUAL
# <DNAStringSet> <DNAStringSetList> <numeric>
#   ML755977.1:177460_ACCCCCCCCCCC/A        ACCCCCCCCCCC                  A   2959.36
# ML755977.1:237488_TAAAAAAAAAAAAA/T      TAAAAAAAAAAAAA                  T   3912.36
# ML755977.1:351610_G/GT                   G                 GT    946.34
# ML755977.1:601795_T/TC                   T                 TC   2291.36
# ML755978.1:83_GAAAAAAAAAAAAAAAAAA/G GAAAAAAAAAAAAAAAAAA                  G  12414.36
# ML755978.1:66627_AGGGGGGG/A            AGGGGGGG                  A   1893.18
# FILTER
# <character>
#   ML755977.1:177460_ACCCCCCCCCCC/A   SORfilter
# ML755977.1:237488_TAAAAAAAAAAAAA/T        PASS
# ML755977.1:351610_G/GT        PASS
# ML755977.1:601795_T/TC   SORfilter
# ML755978.1:83_GAAAAAAAAAAAAAAAAAA/G   SORfilter
# ML755978.1:66627_AGGGGGGG/A   SORfilter

head(info(readvcf.plasmid))
# DataFrame with 6 rows and 18 columns
# AC            AF        AN BaseQRankSum        DP
# <IntegerList> <NumericList> <integer>    <numeric> <integer>
#   ML755977.1:177460_ACCCCCCCCCCC/A                3             1         3        1.280       153
# ML755977.1:237488_TAAAAAAAAAAAAA/T              3             1         3       -1.503       594
# ML755977.1:351610_G/GT                          3             1         3       -0.214       548
# ML755977.1:601795_T/TC                          3             1         3       -0.504       345
# ML755978.1:83_GAAAAAAAAAAAAAAAAAA/G             3             1         3        2.350      1024
# ML755978.1:66627_AGGGGGGG/A                     3             1         3       -0.129       212
# DS       END ExcessHet        FS InbreedingCoeff
# <logical> <integer> <numeric> <numeric>       <numeric>
#   ML755977.1:177460_ACCCCCCCCCCC/A        FALSE        NA        NA     0.000              NA
# ML755977.1:237488_TAAAAAAAAAAAAA/T      FALSE        NA        NA    34.798              NA
# ML755977.1:351610_G/GT                  FALSE        NA        NA     1.249              NA
# ML755977.1:601795_T/TC                  FALSE        NA        NA    18.288              NA
# ML755978.1:83_GAAAAAAAAAAAAAAAAAA/G     FALSE        NA        NA    31.489              NA
# ML755978.1:66627_AGGGGGGG/A             FALSE        NA        NA     0.000              NA
# MLEAC         MLEAF        MQ MQRankSum        QD
# <IntegerList> <NumericList> <numeric> <numeric> <numeric>
#   ML755977.1:177460_ACCCCCCCCCCC/A                3             1     59.86     0.000     25.36
# ML755977.1:237488_TAAAAAAAAAAAAA/T              3             1     59.35     0.000     15.40
# ML755977.1:351610_G/GT                          3             1     60.00     0.000     14.79
# ML755977.1:601795_T/TC                          3             1     59.63     0.000     12.59
# ML755978.1:83_GAAAAAAAAAAAAAAAAAA/G             3             1     49.13    -0.572     28.73
# ML755978.1:66627_AGGGGGGG/A                     2         0.667     60.00     0.000     30.97
# RAW_MQandDP ReadPosRankSum       SOR
# <IntegerList>      <numeric> <numeric>
#   ML755977.1:177460_ACCCCCCCCCCC/A            NA,NA             NA     5.670
# ML755977.1:237488_TAAAAAAAAAAAAA/T          NA,NA         -0.885     0.108
# ML755977.1:351610_G/GT                      NA,NA         -0.041     0.506
# ML755977.1:601795_T/TC                      NA,NA         -0.636     4.144
# ML755978.1:83_GAAAAAAAAAAAAAAAAAA/G         NA,NA             NA     7.321
# ML755978.1:66627_AGGGGGGG/A                 NA,NA          0.979     5.838

write.table(geno(readvcf.plasmid)$DP, "results/readdepth-plasmid.txt", sep = "\t")
write.table(geno(readvcf.strains)$DP, "results/readdepth-strains.txt", sep = "\t")

write.table(geno(readvcf.plasmid)$GT, "results/mutations-plasmid.txt", sep = "\t")
write.table(geno(readvcf.strains)$GT, "results/mutations-strains.txt", sep = "\t")
## here #### 

# Working with genome annotations #### 

sig.hits.plasmid <- rowRanges(readvcf.plasmid)
sig.hits.strains <- rowRanges(readvcf.strains)

my_ld <- 1000

search.plasmid <- flank(sig.hits.plasmid, my_ld, both = TRUE)
search.strains <- flank(sig.hits.strains, my_ld, both = TRUE)

hits.plasmid <- findOverlaps(search.plasmid, gtf.gene.plasmid)
hits.plasmid

hits.plasmid2 <- data.frame(SNP = names(search.plasmid)[queryHits(hits.plasmid)],
                       Gene = gtf.gene.plasmid$ID[subjectHits(hits.plasmid)])

hits.plasmid3 <- findOverlaps(sig.hits.plasmid, gtf.gene.plasmid)

hits.plasmid4 <- data.frame(SNP = names(search.plasmid)[queryHits(hits.plasmid3)],
                            Gene = gtf.gene.plasmid$ID[subjectHits(hits.plasmid3)])
hits.plasmid4

hits.strains <- findOverlaps(search.strains, gtf.gene.strains)
hits.strains

hits.strains2 <- data.frame(SNP = names(search.strains)[queryHits(hits.strains)],
                            Gene = gtf.gene.strains$ID[subjectHits(hits.strains)])

hits.strains3 <- findOverlaps(sig.hits.strains, gtf.gene.strains)

hits.strains4 <- data.frame(SNP = names(search.strains)[queryHits(hits.strains3)],
                            Gene = gtf.gene.strains$ID[subjectHits(hits.strains3)])
hits.strains4

write.table(hits.plasmid4, "results/snpsingene-annotated-plasmid.txt", sep = "\t")
write.table(hits.strains4, "results/snpsingene-annotated-strains.txt", sep = "\t")


save.image("02_Dataanalysis.Rdata")

# Identifying functional consequences of SNPs #### 

snps2search.plasmid <- subsetByOverlaps(sig.hits.plasmid, gtf.gene.plasmid)
snps2search.strains <- subsetByOverlaps(sig.hits.strains, gtf.gene.strains)

# genome.plasmid <- FaFile("data/pINT03-XYL123.fa")
txdb.strains <- makeTxDbFromGFF("data/GCA_009372015.1_YarliW29_genomic.gff",
                           format = "gff", organism = "Yarrowia lipolytica",
                           dataSource = "NCBI")

txdb.plasmid <- makeTxDbFromGFF("data/pINT03-XYL123.gff",
                                format = "gff", organism = "Yarrowia lipolytica",
                                dataSource = "NCBI")


genes(txdb.strains)

exonsBy(txdb.strains, by = "tx")

# error here. 

# tomorrow - output the list of 20 and work with that. 
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

save.image("Mar04-1618.Rdata")



# ISSUE Code Start #### 

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

