# read depth analysis 


# setup libraries #### 
library(magrittr)
library(ggplot2)
library(ggpubr)

rm(list = ls())

# data import #### 
getwd()
setwd("D:/Box Sync/02 - RNA seq work/2101-YL-Evolution (Sangdo)/R analysis")

depth.xev <- read.table("data/coverage_depth/XEV_GTATGTTC-TTCCTGTT_aligned.sorted.depth.txt")
depth.xev <- depth.xev[,2:3]
colnames(depth.xev) <- c("Pos", "XEV_Depth")
head(depth.xev)

depth.gev <- read.table("data/coverage_depth/GEV_ACGCACCT-CCTTCACC_aligned.sorted.depth.txt")
depth.gev <- depth.gev[,2:3]
colnames(depth.gev) <- c("Pos", "GEV_Depth")
head(depth.gev)

depth.x12 <- read.table("data/coverage_depth/X123_CGCTATGT-GTGTCGGA_aligned.sorted.depth.txt")
depth.x12 <- depth.x12[,2:3]
colnames(depth.x12) <- c("Pos", "X12_Depth")
head(depth.x12)

plasmid.features <- read.delim("data/plasmid_features.txt")

# Combine all three samples #### 
'And remove any errors and NA values '
summary(depth.xev$Pos)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1    3307    6612    6612    9918   13224 
summary(depth.gev$Pos)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1    3315    6618    6617    9921   13224 
summary(depth.x12$Pos)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 181    3449    6708    6707    9966   13224 

summary(depth.xev$XEV_Depth)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0    3951    5428    4694    5779    6560 
summary(depth.gev$GEV_Depth)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0    1795    2063    1971    2421    2967 
summary(depth.x12$X12_Depth)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0   665.0   788.0   723.7   843.0   995.0 

depth.all.1 <- merge.data.frame(depth.gev, depth.xev, by = "Pos", all = TRUE, sort = TRUE, no.dups = TRUE)
depth.all <- merge.data.frame(depth.all.1, depth.x12, by = "Pos", all = TRUE, sort = TRUE, no.dups = TRUE)

rm(depth.all.1)

summary(depth.all$Pos)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 1    3307    6612    6612    9918   13224 
summary(depth.all$XEV_Depth)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0    3951    5428    4694    5779    6560 
summary(depth.all$GEV_Depth)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NAs 
#       0    1795    2063    1971    2421    2967      11 
summary(depth.all$X12_Depth)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NAs 
#     0.0   665.0   788.0   723.7   843.0   995.0     190

# fixing NA's
depth.all[is.na(depth.all)] = 0

summary(depth.all$GEV_Depth)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max
#      0    1794    2062    1970    2421    2967 
summary(depth.all$X12_Depth)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max
#    0.0   656.0   787.0   713.3   842.0   995.0

# Data cleanup of low read depths #### 

jpeg(filename = "results/readdepth_distribution.jpeg", width = 10, 
     height = 10, units = "in", res = 600)
# x11(w = 10, h = 10)
a <- ggplot(depth.all, aes(x = X12_Depth)) + geom_density() +
  ggtitle("Density distribution - X123") + xlab("Read Depth") +
  theme(plot.title = element_text(size=12, face="bold", 
                                  margin = margin(20, 20, 20, 20))) 
b <- ggplot(depth.all, aes(x = XEV_Depth)) + geom_density() +
  ggtitle("Density distribution - XEV") + xlab("Read Depth") +
  theme(plot.title = element_text(size=12, face="bold", 
                                  margin = margin(20, 20, 20, 20))) 
c <- ggplot(depth.all, aes(x = GEV_Depth)) + geom_density() +
  ggtitle("Density distribution - GEV") + xlab("Read Depth") +
  theme(plot.title = element_text(size=12, face="bold", 
                                  margin = margin(20, 20, 20, 20))) 
ggarrange(a, b, c, ncol = 1, nrow = 3)
dev.off()

sum(depth.all$X12_Depth > 50)/length(depth.all$X12_Depth)
# 0.9766334
sum(depth.all$XEV_Depth > 100)/length(depth.all$XEV_Depth)
# 0.9726255
sum(depth.all$GEV_Depth > 100)/length(depth.all$GEV_Depth)
# 0.97134

depth.all$X12_Depth[depth.all$X12_Depth <= 50] <- 0 
depth.all$XEV_Depth[depth.all$XEV_Depth <= 100] <- 0 
depth.all$GEV_Depth[depth.all$GEV_Depth <= 100] <- 0 


# Read depth ratios #### 

depth.all$XEV.X12 = depth.all$XEV_Depth/depth.all$X12_Depth
depth.all$GEV.X12 = depth.all$GEV_Depth/depth.all$X12_Depth
depth.all$XEV.GEV = depth.all$XEV_Depth/depth.all$GEV_Depth

depth.all[is.na(depth.all)] = 0

# Add plasmid features for plotting #### 

plasmid.features
# Feature Start   End
# 1             Lip7   186  1292
# 2            GPM1p  3235  4078
# 3             XYL1  4079  5035
# 4            GPM1t  5036  5335
# 5            TDH1p  5336  6302
# 6             XYL2  6303  7394
# 7            TDH1t  7395  7694
# 8            FBA1p  7695  8526
# 9             XYL3  8527 10398
# 10           FBA1t 10399 10698
# 11  cyc terminator 10711 10997
# 12 ColE1 origin(1) 11365 12047
# 13            AmpR 12145 12804

depth.all$Features <- c("Others")
depth.all$Features[186:1292] <- c("Lip7")
depth.all$Features[3235:4078] <- c("GPM1p")
depth.all$Features[4079:5035] <- c("XYL1")
depth.all$Features[5036:5335] <- c("GPM1t")
depth.all$Features[5336:6302] <- c("TDH1p")
depth.all$Features[6303:7394] <- c("XYL2")
depth.all$Features[7395:7694] <- c("TDH1t")
depth.all$Features[7695:8526] <- c("FBA1p")
depth.all$Features[8527:10398] <- c("XYL3")
depth.all$Features[10399:10698] <- c("FBA1t")
depth.all$Features[10711:10997] <- c("cyc_term")
depth.all$Features[12145:12804] <- c("AmpR")

depth.all$Gene <- c("Others")
depth.all$Gene[4079:5035] <- c("XYL1")
depth.all$Gene[6303:7394] <- c("XYL2")
depth.all$Gene[8527:10398] <- c("XYL3")


# Plotting ####

# x11(w = 15, h = 10)
jpeg(filename = "results/readdepth_ratios.jpeg", width = 15,
     height = 10, units = "in", res = 600)
a <- ggplot(depth.all, aes(x = Pos, y = XEV.X12, color = Gene)) + 
  geom_line() + ylim(0,10) +
  scale_color_manual(values=c("#999999", "#4682B4", "#B4464B", "#B4AF46")) + 
  xlab("Plasmid coordinates") + ylab("Ratio") + ggtitle("Relative read depth: XEV / X123") +
  theme(plot.title = element_text(size=12, face="bold", margin = margin(20, 20, 20, 20))) + 
  guides(colour = guide_legend(override.aes = list(size=4))) +
  geom_vline(xintercept = 4079, linetype="dashed", color = "#4682B4", size=0.75) + 
  geom_vline(xintercept = 5035, linetype="dashed", color = "#4682B4", size=0.75) + 
  geom_vline(xintercept = 6303, linetype="dashed", color = "#B4464B", size=0.75) + 
  geom_vline(xintercept = 7394, linetype="dashed", color = "#B4464B", size=0.75) + 
  geom_vline(xintercept = 8527, linetype="dashed", color = "#B4AF46", size=0.75) + 
  geom_vline(xintercept = 10398, linetype="dashed", color = "#B4AF46", size=0.75)

b <- ggplot(depth.all, aes(x = Pos, y = GEV.X12, color = Gene)) + 
  geom_line() + ylim(0,10) +
  scale_color_manual(values=c("#999999", "#4682B4", "#B4464B", "#B4AF46")) +  
  xlab("Plasmid coordinates") + ylab("Ratio") + ggtitle("Relative read depth: GEV / X123") +
  theme(plot.title = element_text(size=12, face="bold", margin = margin(20, 20, 20, 20))) + 
  guides(colour = guide_legend(override.aes = list(size=4))) +
  geom_vline(xintercept = 4079, linetype="dotted", color = "#4682B4", size=0.75) + 
  geom_vline(xintercept = 5035, linetype="dotted", color = "#4682B4", size=0.75) + 
  geom_vline(xintercept = 6303, linetype="dotted", color = "#B4464B", size=0.75) + 
  geom_vline(xintercept = 7394, linetype="dotted", color = "#B4464B", size=0.75) + 
  geom_vline(xintercept = 8527, linetype="dotted", color = "#B4AF46", size=0.75) + 
  geom_vline(xintercept = 10398, linetype="dotted", color = "#B4AF46", size=0.75)

c <- ggplot(depth.all, aes(x = Pos, y = XEV.GEV, color = Gene)) + 
  geom_line() + ylim(0,10) +
  scale_color_manual(values=c("#999999", "#4682B4", "#B4464B", "#B4AF46")) + 
  xlab("Plasmid coordinates") + ylab("Ratio") + ggtitle("Relative read depth: XEV / GEV") +
  theme(plot.title = element_text(size=12, face="bold", margin = margin(20, 20, 20, 20))) + 
  guides(colour = guide_legend(override.aes = list(size=4))) +
  geom_vline(xintercept = 4079, linetype="dashed", color = "#4682B4", size=0.75) + 
  geom_vline(xintercept = 5035, linetype="dashed", color = "#4682B4", size=0.75) + 
  geom_vline(xintercept = 6303, linetype="dashed", color = "#B4464B", size=0.75) + 
  geom_vline(xintercept = 7394, linetype="dashed", color = "#B4464B", size=0.75) + 
  geom_vline(xintercept = 8527, linetype="dashed", color = "#B4AF46", size=0.75) + 
  geom_vline(xintercept = 10398, linetype="dashed", color = "#B4AF46", size=0.75)

ggarrange(a, b, c, ncol = 1, nrow = 3)
dev.off()

jpeg(filename = "results/readdepth_ratio_distribution.jpeg", width = 10,
     height = 10, units = "in", res = 600)
# x11(w = 10, h = 10)
a <- ggplot(depth.all, aes(x = XEV.X12, color = Gene)) + geom_density() + 
  xlab("Relative Read Depth") + xlim(0,10) + ggtitle("Read Depth Distribution: XEV / X123") +
  theme(plot.title = element_text(size=12, face="bold", margin = margin(20, 20, 20, 20))) 
b <- ggplot(depth.all, aes(x = GEV.X12, color = Gene)) + geom_density() +
  xlab("Relative Read Depth") + xlim(0,10) + ggtitle("Read Depth Distribution: GEV / X123") +
  theme(plot.title = element_text(size=12, face="bold", margin = margin(20, 20, 20, 20))) 
c <- ggplot(depth.all, aes(x = XEV.GEV, color = Gene)) + geom_density() +
  xlab("Relative Read Depth") + xlim(0,10) + ggtitle("Read Depth Distribution: XEV / GEV") +
  theme(plot.title = element_text(size=12, face="bold", margin = margin(20, 20, 20, 20))) 

ggarrange(a, b, c, ncol = 1, nrow = 3)
dev.off()

# Read depth calculation #### 

xyl1 <- depth.all[depth.all$Features == "XYL1",]
xyl2 <- depth.all[depth.all$Features == "XYL2",]
xyl3 <- depth.all[depth.all$Features == "XYL3",]

summary(xyl1)
# Pos         GEV_Depth      XEV_Depth      X12_Depth        XEV.X12         GEV.X12         XEV.GEV        Features             Gene          
# Min.   :4079   Min.   :1765   Min.   :4078   Min.   :545.0   Min.   :6.353   Min.   :2.670   Min.   :2.178   Length:957         Length:957        
# 1st Qu.:4318   1st Qu.:2420   1st Qu.:5562   1st Qu.:768.0   1st Qu.:6.996   1st Qu.:3.088   1st Qu.:2.257   Class :character   Class :character  
# Median :4557   Median :2482   Median :5699   Median :789.0   Median :7.216   Median :3.162   Median :2.294   Mode  :character   Mode  :character  
# Mean   :4557   Mean   :2474   Mean   :5678   Mean   :792.1   Mean   :7.185   Mean   :3.131   Mean   :2.296                                        
# 3rd Qu.:4796   3rd Qu.:2535   3rd Qu.:5826   3rd Qu.:819.0   3rd Qu.:7.419   3rd Qu.:3.222   3rd Qu.:2.329                                        
# Max.   :5035   Max.   :2676   Max.   :6008   Max.   :909.0   Max.   :8.485   Max.   :3.687   Max.   :2.440                                     
# 
summary(xyl2)
# Pos         GEV_Depth      XEV_Depth      X12_Depth        XEV.X12         GEV.X12         XEV.GEV        Features             Gene          
# Min.   :6303   Min.   :1949   Min.   :4562   Min.   :679.0   Min.   :6.190   Min.   :2.658   Min.   :2.174   Length:1092        Length:1092       
# 1st Qu.:6576   1st Qu.:2446   1st Qu.:5747   1st Qu.:828.0   1st Qu.:6.611   1st Qu.:2.832   1st Qu.:2.277   Class :character   Class :character  
# Median :6848   Median :2535   Median :5870   Median :868.5   Median :6.783   Median :2.913   Median :2.315   Mode  :character   Mode  :character  
# Mean   :6848   Mean   :2540   Mean   :5880   Mean   :866.9   Mean   :6.794   Mean   :2.933   Mean   :2.318                                        
# 3rd Qu.:7121   3rd Qu.:2621   3rd Qu.:6000   3rd Qu.:908.0   3rd Qu.:6.946   3rd Qu.:3.033   3rd Qu.:2.362                                        
# Max.   :7394   Max.   :2897   Max.   :6517   Max.   :972.0   Max.   :7.536   Max.   :3.272   Max.   :2.521
# 
summary(xyl3)
# Pos          GEV_Depth      XEV_Depth      X12_Depth        XEV.X12         GEV.X12         XEV.GEV        Features             Gene          
# Min.   : 8527   Min.   :1772   Min.   :4370   Min.   :643.0   Min.   :6.168   Min.   :2.583   Min.   :2.178   Length:1872        Length:1872       
# 1st Qu.: 8995   1st Qu.:2341   1st Qu.:5591   1st Qu.:813.0   1st Qu.:6.680   1st Qu.:2.776   1st Qu.:2.314   Class :character   Class :character  
# Median : 9462   Median :2407   Median :5716   Median :835.0   Median :6.830   Median :2.896   Median :2.361   Mode  :character   Mode  :character  
# Mean   : 9462   Mean   :2406   Mean   :5694   Mean   :835.6   Mean   :6.820   Mean   :2.882   Mean   :2.369                                        
# 3rd Qu.: 9930   3rd Qu.:2490   3rd Qu.:5832   3rd Qu.:860.2   3rd Qu.:6.970   3rd Qu.:2.983   3rd Qu.:2.427                                        
# Max.   :10398   Max.   :2664   Max.   :6116   Max.   :934.0   Max.   :8.094   Max.   :3.389   Max.   :2.618 
# 
# 
summary(depth.all)
# Pos          GEV_Depth      XEV_Depth      X12_Depth        XEV.X12         GEV.X12         XEV.GEV        Features             Gene          
# Min.   :    1   Min.   :   0   Min.   :   0   Min.   :  0.0   Min.   :0.000   Min.   :0.000   Min.   :0.000   Length:13224       Length:13224      
# 1st Qu.: 3307   1st Qu.:1794   1st Qu.:3951   1st Qu.:656.0   1st Qu.:6.215   1st Qu.:2.478   1st Qu.:2.190   Class :character   Class :character  
# Median : 6612   Median :2062   Median :5428   Median :787.0   Median :6.738   Median :2.773   Median :2.310   Mode  :character   Mode  :character  
# Mean   : 6612   Mean   :1969   Mean   :4692   Mean   :713.1   Mean   :  Inf   Mean   :  Inf   Mean   :  Inf                                        
# 3rd Qu.: 9918   3rd Qu.:2421   3rd Qu.:5779   3rd Qu.:842.0   3rd Qu.:7.029   3rd Qu.:3.025   3rd Qu.:2.468                                        
# Max.   :13224   Max.   :2967   Max.   :6560   Max.   :995.0   Max.   :  Inf   Max.   :  Inf   Max.   :  Inf   

write.table(summary(xyl1), "results/read_depth_xyl1.txt")
write.table(summary(xyl2), "results/read_depth_xyl2.txt")
write.table(summary(xyl3), "results/read_depth_xyl3.txt")
write.table(summary(depth.all), "results/read_depth_all.txt")

save.image("202105_Readdepth_analysis.Rdata")

sessionInfo()
# R version 4.0.5 (2021-03-31)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 19042)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
# [5] LC_TIME=English_United States.1252    
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] ggpubr_0.4.0                gplots_3.1.1                scater_1.18.6               SingleCellExperiment_1.12.0 ggplot2_3.3.3              
# [6] snpStats_1.40.0             Matrix_1.3-2                survival_3.2-10             VariantAnnotation_1.36.0    magrittr_2.0.1             
# [11] SummarizedExperiment_1.20.0 MatrixGenerics_1.2.1        matrixStats_0.58.0          Rsamtools_2.6.0             Biostrings_2.58.0          
# [16] XVector_0.30.0              GenomicFeatures_1.42.2      AnnotationDbi_1.52.0        Biobase_2.50.0              GenomicRanges_1.42.0       
# [21] GenomeInfoDb_1.26.4         IRanges_2.24.1              S4Vectors_0.28.1            BiocGenerics_0.36.0        
# 
# loaded via a namespace (and not attached):
#   [1] ggbeeswarm_0.6.0          colorspace_2.0-0          ggsignif_0.6.1            ellipsis_0.3.1            rio_0.5.26                scuttle_1.0.4            
# [7] BiocNeighbors_1.8.2       rstudioapi_0.13           farver_2.1.0              bit64_4.0.5               fansi_0.4.2               xml2_1.3.2               
# [13] splines_4.0.5             sparseMatrixStats_1.2.1   cachem_1.0.4              broom_0.7.6               dbplyr_2.1.0              BiocManager_1.30.12      
# [19] compiler_4.0.5            httr_1.4.2                backports_1.2.1           assertthat_0.2.1          fastmap_1.1.0             BiocSingular_1.6.0       
# [25] prettyunits_1.1.1         tools_4.0.5               rsvd_1.0.3                gtable_0.3.0              glue_1.4.2                GenomeInfoDbData_1.2.4   
# [31] dplyr_1.0.5               rappdirs_0.3.3            tinytex_0.31              Rcpp_1.0.6                carData_3.0-4             cellranger_1.1.0         
# [37] vctrs_0.3.6               rtracklayer_1.49.5        DelayedMatrixStats_1.12.3 xfun_0.22                 stringr_1.4.0             openxlsx_4.2.3           
# [43] beachmat_2.6.4            lifecycle_1.0.0           irlba_2.3.3               gtools_3.8.2              rstatix_0.7.0             XML_3.99-0.6             
# [49] zlibbioc_1.36.0           scales_1.1.1              BSgenome_1.58.0           hms_1.0.0                 curl_4.3                  memoise_2.0.0            
# [55] gridExtra_2.3             biomaRt_2.46.3            stringi_1.5.3             RSQLite_2.2.5             caTools_1.18.2            zip_2.1.1                
# [61] BiocParallel_1.24.1       rlang_0.4.10              pkgconfig_2.0.3           bitops_1.0-6              lattice_0.20-41           purrr_0.3.4              
# [67] GenomicAlignments_1.26.0  labeling_0.4.2            cowplot_1.1.1             bit_4.0.4                 tidyselect_1.1.0          R6_2.5.0                 
# [73] generics_0.1.0            DelayedArray_0.16.3       DBI_1.1.1                 pillar_1.5.1              haven_2.4.1               foreign_0.8-81           
# [79] withr_2.4.1               abind_1.4-5               RCurl_1.98-1.3            tibble_3.1.0              crayon_1.4.1              car_3.0-10               
# [85] KernSmooth_2.23-18        utf8_1.2.1                BiocFileCache_1.14.0      viridis_0.5.1             progress_1.2.2            readxl_1.3.1             
# [91] grid_4.0.5                data.table_1.14.0         blob_1.2.1                forcats_0.5.1             digest_0.6.27             tidyr_1.1.3              
# [97] openssl_1.4.3             munsell_0.5.0             beeswarm_0.3.1            viridisLite_0.3.0         vipor_0.4.5               askpass_1.1 
