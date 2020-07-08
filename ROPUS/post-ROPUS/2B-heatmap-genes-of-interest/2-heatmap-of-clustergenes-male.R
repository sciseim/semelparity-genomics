setwd("~/Dropbox/seimlab at NNU/antechinus genome project/--data and analysis--/---late_2019---/ROPUS_data/flavipes_RNAseq/FCM/") # OSX
owd <- getwd()

# load all normalised counts
counts.LOG <- readRDS("../normalized_counts.rds")

# the.tissue <- "kidney"
# the.tissue <- "liver"
# the.tissue <- "skeletalmuscle"
# the.tissue <- "cerebrum"
# the.tissue <- "testis"

the.sex <- "male"



targetgene <- c("NEMP2","SH3TC1","VCPKMT","POLR1D","NEUROD6","PPP1R3G","BRK1","TMEM271","TRAPPC6A","PANX1","TEFM","ATP1B4","APOB","C3orf80","C12orf43","TWISTNB","STRA8","TMC3","ANKS4B","VN1R2","COMMD5","TRAPPC2L","FMC1","EPHA1","CCDC167","FAM207A","SEPSECS","YRDC","FANCM","IFIT5","MGME1","GPR34","SMIM26","CETP","H2AFV","SMCO1","TLR9","TRAF3IP3","ENPP1","ADGRV1","RRNAD1","ATAD5","ZP2","BAX","TMEM256","TMEM37","C6orf58","TMEM14A","FKBP5","TANGO6","LYPLAL1","MGST3","KHK","SELENOH","CCDC32","EDAR","HBQ1","CAMK2N1","AOPEP","PPIH","PIGN","PEMT","PROX2","DNASE1","LSM3","CYB561D2","ACP6","VSIG2","TRIP13","JAKMIP3","WNT4","NT5DC1","COX14","POLR3K","MRI1","COA5","NDUFC1","AC005697.1","SAC3D1","SNRPE","MOCS3","UPK3B","KCNG4","PMVK","CHCHD1","ZNF444","MRPS18C","NIPA2","PEX2","COMMD1","HTD2","MRPL46","LSM5","METTL6","TSPO","TMEM126A","AP4S1","TARBP2","CNIH4","DCPS","NDUFAF5","RABEPK","PLGRKT","TMEM258","MDP1","CNIH1","HGFAC","SIVA1","CXorf40A","TMEM251","AC084337.2","KCNH5","SDCBP2","PTGS2","POP7","APOM","CNTD1","CSTB","PSPH","IDNK","CMC1","HORMAD1","OR51B6","PMM2","HECA","TRMT6","KLHL28","ZBTB8A","TMEM18","PTCHD1","ZNF367","CYSLTR1","ITGA11","CCR9","HECTD2","EIF2B3","GNG12","MRPL41","CCDC85B","ABT1","MBLAC1","LYPD6","HDDC2","NPFF","WNT10A","TRIML1","PRICKLE1","CTHRC1","FBXL4","BET1","TMEM9B","MRPS10","COX11","SLC17A5","DNAJC30","MORN3","SSC4D","BMT2","FDX1","PAFAH2","LZTFL1","KCNB2","SPATA5L1","HEPH","PSMG4","ODR4","CPA2","RALGPS1","SERP1","MRPS36","CBY1","SPC25","GINS4","NDUFAF8","NTHL1","DNLZ","RPIA","NFE2","EBP","UQCC3","TREH","DCXR","TRMT44","VSX1","BLOC1S6","COMTD1","BCO1","FSD2","SEC61G","AGRP","TSPAN15","TIMM8B","ARSE","MSRB2","PFDN4")

# OXPHOS genes (from http://software.broadinstitute.org/gsea/msigdb/cards/HALLMARK_OXIDATIVE_PHOSPHORYLATION)
targetgene <- c("ABCB7","ACAA1","ACAA2","ACADM","ACADSB","ACADVL","ACAT1","ACO2","AFG3L2","AIFM1","ALAS1","ALDH6A1","ATP1B1","ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E","ATP5PB","ATP5MC1","ATP5MC2","ATP5MC3","ATP5PD","ATP5ME","ATP5PF","ATP5MF","ATP5MG","ATP5PO","ATP6AP1","ATP6V0B","ATP6V0C","ATP6V0E1","ATP6V1C1","ATP6V1D","ATP6V1E1","ATP6V1F","ATP6V1G1","ATP6V1H","BAX","BCKDHA","BDH2","MPC1","CASP7","COX10","COX11","COX15","COX17","COX4I1","COX5A","COX5B","COX6A1","COX6B1","COX6C","COX7A2","COX7A2L","COX7B","COX7C","COX8A","CPT1A","CS","CYB5A","CYB5R3","CYC1","CYCS","DECR1","DLAT","DLD","DLST","ECH1","ECHS1","ECI1","ETFA","ETFB","ETFDH","FDX1","FH","FXN","GLUD1","GOT2","GPI","GPX4","GRPEL1","HADHA","HADHB","HCCS","HSD17B10","HSPA9","HTRA2","IDH1","IDH2","IDH3A","IDH3B","IDH3G","IMMT","ISCA1","ISCU","LDHA","LDHB","LRPPRC","MAOB","MDH1","MDH2","MFN2","MGST3","MRPL11","MRPL15","MRPL34","MRPL35","MRPS11","MRPS12","MRPS15","MRPS22","MRPS30","MTRF1","MTRR","MTX2","NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8","NDUFA9","NDUFAB1","NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFC1","NDUFC2","NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS6","NDUFS7","NDUFS8","NDUFV1","NDUFV2","NNT","NQO2","OAT","OGDH","OPA1","OXA1L","PDHA1","PDHB","PDHX","PDK4","PDP1","PHB2","PHYH","PMPCA","POLR2F","POR","PRDX3","RETSAT","RHOT1","RHOT2","SDHA","SDHB","SDHC","SDHD","SLC25A11","SLC25A12","SLC25A20","SLC25A3","SLC25A4","SLC25A5","SLC25A6","SUCLA2","SUCLG1","SUPV3L1","SURF1","TCIRG1","TIMM10","TIMM13","TIMM17A","TIMM50","TIMM8B","TIMM9","TOMM22","TOMM70","UQCR10","UQCR11","UQCRB","UQCRC1","UQCRC2","UQCRFS1","UQCRH","UQCRQ","VDAC1","VDAC2","VDAC3")

# keep only samples samples used in FCM
# male
counts.by.samling.date.male <- counts.LOG[,c("MN21818Li","MB22918Li","MB22918Li2","MB25918Li","MB25918Li2","MB021018Li","MB021018Li2","MB041018Li",
                                             "MN21818Sm","MB22918Sm","MB22918Sm2","MB25918Sm","MB25918Sm2","MB021018Sm","MB021018Sm2","MB041018Sm",
                                             "MN21818Ki","MB22918Ki","MB22918Ki2","MB25918Ki","MB25918Ki2","MB021018Ki","MB021018Ki2","MB041018Ki",
                                             "MN21818Ce","MB22918Ce","MB22918Ce2","MB25918Ce2","MB021018Ce","MB021018Ce2","MB041018Ce",                                        "MN21818Te","MB22918Te","MB22918Te2","MB25918Te","MB25918Te2","MB021018Te","MB021018Te2","MB041018Te")]


# or a subset 
# LIVER
if(the.tissue == "liver")
{ counts.by.samling.date.male <- counts.LOG[,c("MN21818Li","MB22918Li","MB22918Li2","MB25918Li","MB25918Li2","MB021018Li","MB021018Li2","MB041018Li")] }
# SKELETAL MUSCLE
if(the.tissue == "skeletalmuscle")
{ counts.by.samling.date.male <- counts.LOG[,c("MN21818Sm","MB22918Sm","MB22918Sm2","MB25918Sm","MB25918Sm2","MB021018Sm","MB021018Sm2","MB041018Sm")] }
# KIDNEY
if(the.tissue == "kidney")
{ counts.by.samling.date.male <- counts.LOG[,c("MN21818Ki","MB22918Ki","MB22918Ki2","MB25918Ki","MB25918Ki2","MB021018Ki","MB021018Ki2","MB041018Ki")] }
# CEREBRUM
if(the.tissue == "cerebrum")
{ counts.by.samling.date.male <- counts.LOG[,c("MN21818Ce","MB22918Ce","MB22918Ce2","MB25918Ce2","MB021018Ce","MB021018Ce2")] }
# TESTIS
if(the.tissue == "testis")
{ counts.by.samling.date.male <- counts.LOG[,c("MN21818Te","MB22918Te","MB22918Te2","MB25918Te","MB25918Te2","MB021018Te","MB021018Te2","MB041018Te")] }








# assign to temp DF
counts.by.samling.date.TEMP <- counts.by.samling.date.male


# subset by genes of interest
# counts.df2 <- 2^(  ( counts.by.samling.date.TEMP[which(row.names(counts.by.samling.date.TEMP) %in% targetgene),]) )-0.25
# class(counts.df2) # df
counts.df2 <-   counts.by.samling.date.TEMP[which(row.names(counts.by.samling.date.TEMP) %in% targetgene),] 
class(counts.df2) # df

# sort by input order
# library(dplyr)
# counts.df2 <- counts.df2[targetgene,]

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# HEATMAP
library("pheatmap")
mat <- as.matrix(counts.df2)



# ######################################################################################
# HEAT MAPS
# ######################################################################################
#Load latest version of heatmap.3 function
# source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
source("heatmap.3.R")
#
#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}

# heatmap colour parameters
cols <- colorRampPalette(c("green", "black", "red"))(n = 256) # classical
cols <- rev(colorRampPalette(brewer.pal(11, "RdBu")) (256))

# documentation: https://rdrr.io/github/lwaldron/LeviRmisc/man/heatmap.3.html
# pdf("TEST.pdf",height=10)
# h <- heatmap.3((mat), col=cols, scale="row", trace="none",density.info="none", dendrogram="none",Colv=FALSE, Rowv=TRUE, key=TRUE,cexRow=0.5,breaks=seq(-3,3,6/256)   )
# print(h)
# dev.off()


# targetgene <- c("BPGM","FBP2","GPI","LDHA","PFKM","PGK2","PKM","GPD1","GPD2","SLC2A5","SLC37A4","PRKAG3","STK11","SOCS1","SOCS2")

counts.df2["PKM",]
counts.df2["BPGM",]

# ColSideColors=column_annotation,RowSideColors=rlab

the.file.name <- paste(the.sex,"-",the.tissue,"_heatmap.pdf",sep="")
# pdf(the.file.name,height=30)

# pdf(the.file.name)
# h <- heatmap.3((counts.df2), col=cols, scale="row", trace="none",density.info="none", dendrogram="none",Colv=FALSE, Rowv=TRUE, key=TRUE,cexRow=0.5,breaks=seq(-3,3,6/256)   )
# print(h)
# dev.off()

??heatmap.3
pdf(the.file.name,height=20)
heatmap.3((counts.df2), col=cols, scale="row", trace="none",density.info="none",denscol="black", dendrogram="none",Colv=FALSE, Rowv=FALSE, key=TRUE,cexRow=0.5,breaks=seq(-3,3,6/256)   )
dev.off()