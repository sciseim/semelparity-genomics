rm(list=ls()) # reset R  
setwd("~/Dropbox/seimlab at NNU/antechinus genome project/--data and analysis--/---late_2019---/ROPUS_data/flavipes_RNAseq/FCM/") # OSX
owd <- getwd()

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# load all normalised counts
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
counts.LOG <- readRDS("../normalized_counts.rds")


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# SET VARIABLES
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
the.sex <- "female"





# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# KEEP ONLY MALE OR FEMALE SAMPLES
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# rationale: keep only samples samples used in the FCM!


# male
counts.by.samling.date.male <- counts.LOG[,c("MN21818Li","MB22918Li","MB22918Li2","MB25918Li","MB25918Li2","MB021018Li","MB021018Li2","MB041018Li",
                                             "MN21818Sm","MB22918Sm","MB22918Sm2","MB25918Sm","MB25918Sm2","MB021018Sm","MB021018Sm2","MB041018Sm",
                                             "MN21818Ki","MB22918Ki","MB22918Ki2","MB25918Ki","MB25918Ki2","MB021018Ki","MB021018Ki2","MB041018Ki",
                                             "MN21818Ce","MB22918Ce","MB22918Ce2","MB25918Ce2","MB021018Ce","MB021018Ce2","MB041018Ce",                                        "MN21818Te","MB22918Te","MB22918Te2","MB25918Te","MB25918Te2","MB021018Te","MB021018Te2","MB041018Te")]
# or a subset 
#@ counts.by.samling.date.male <- counts.LOG[,c("MN21818Li","MB22918Li","MB22918Li2","MB25918Li","MB25918Li2","MB021018Li","MB021018Li2","MB041018Li")]

counts.by.samling.date.female <- counts.LOG[,c("FN21818Li","FN21818Li2","FB22918Li","FB25918Li","FB021018Li2",
                                               "FN21818Sm","FN21818Sm2","FN21818Sm3","FB22918Sm","FB25918Sm","FB21018Sm2",
                                               "FN21818Ki","FN21818Ki2","FN21818Ki3","FB22918Ki","FB25918Ki","FB21018Ki2",
                                               "FN21818Ce","FN21818Ce2","FB22918Ce","FB25918Ce","FB021018Ce2",
                                               "FN21818Ov","FN21818Ov2","FB22918Ov","FB25918Ov","FB21018Ov2")]



# assign to one DF per tissue per sex
Li.countsM <- counts.by.samling.date.male[,  grep("Li", names(counts.by.samling.date.male) )  ]
Sm.countsM <- counts.by.samling.date.male[,  grep("Sm|Sk", names(counts.by.samling.date.male) )  ]
Ki.countsM <- counts.by.samling.date.male[,  grep("Ki", names(counts.by.samling.date.male) )  ]
Ce.countsM <- counts.by.samling.date.male[,  grep("Ce", names(counts.by.samling.date.male) )  ]
Re.countsM <- counts.by.samling.date.male[,  grep("Te", names(counts.by.samling.date.male) )  ]
Li.countsF <- counts.by.samling.date.female[,  grep("Li", names(counts.by.samling.date.female) )  ]
Sm.countsF <- counts.by.samling.date.female[,  grep("Sm|Sk", names(counts.by.samling.date.female) )  ]
Ki.countsF <- counts.by.samling.date.female[,  grep("Ki", names(counts.by.samling.date.female) )  ]
Ce.countsF <- counts.by.samling.date.female[,  grep("Ce", names(counts.by.samling.date.female) )  ]
Re.countsF <- counts.by.samling.date.female[,  grep("Ov", names(counts.by.samling.date.female) )  ]


# male times
LiMa.time=c("21818","22918","22918","25918","25918","021018","021018","041018")# 10 in our example
SmMa.time=c("21818","22918","22918","25918","25918","021018","021018","041018")
KiMa.time=c("21818","22918","22918","25918","25918","021018","021018","041018")
CeMa.time=c("21818","22918","22918","25918","021018","021018","041018")
ReMa.time=c("21818","22918","22918","25918","25918","021018","021018","041018")
# female times
LiFe.time=c("21818","21818","22918","25918","021018")
SmFe.time=c("21818","21818","21818","22918","25918","021018")
KiFe.time=c("21818","21818","21818","22918","25918","021018")
CeFe.time=c("21818","21818","22918","25918","021018")
ReFe.time=c("21818","21818","22918","25918","021018")
# add to design
LiMa.time=c("other","other","other","other","other","moribund","moribund","moribund") # 10 in our example
SmMa.time=c("other","other","other","other","other","moribund","moribund","moribund")
KiMa.time=c("other","other","other","other","other","moribund","moribund","moribund")
CeMa.time=c("other","other","other","other","moribund","moribund","moribund")
ReMa.time=c("other","other","other","other","other","moribund","moribund","moribund")
# female times
LiFe.time=c("other","other","other","other","moribund")
SmFe.time=c("other","other","other","other","other","moribund")
KiFe.time=c("other","other","other","other","other","moribund")
CeFe.time=c("other","other","other","other","moribund")
ReFe.time=c("other","other","other","other","moribund")


# assign into tissue vectors per sex
if(the.sex == "male"){
  them.tissues <- c("cerebrum","kidney","liver","skeletalmuscle","testis")}
if(the.sex == "female"){
  them.tissues <- c("cerebrum","kidney","liver","ovary","skeletalmuscle")}


for (i in them.tissues) {
  the.tissue <- i
  # START OF BIG LOOP 
  
  if(the.tissue == "liver" & the.sex == "male")
  {
    # testing
    the.timepoint <- LiMa.time
    counts.SUBSET <- Li.countsM
  }
  if(the.tissue == "cerebrum" & the.sex == "male")
  {
    # testing
    the.timepoint <- CeMa.time
    counts.SUBSET <- Ce.countsM
  }
  if(the.tissue == "kidney" & the.sex == "male")
  {
    # testing
    the.timepoint <- KiMa.time
    counts.SUBSET <- Ki.countsM
  }
  if(the.tissue == "testis" & the.sex == "male")
  {
    # testing
    the.timepoint <- ReMa.time
    counts.SUBSET <- Re.countsM
  }
  if(the.tissue == "skeletalmuscle" & the.sex == "male")
  {
    # testing
    the.timepoint <- SmMa.time
    counts.SUBSET <- Sm.countsM
  }
  if(the.tissue == "liver" & the.sex == "female")
  {
    # testing
    the.timepoint <- LiFe.time
    counts.SUBSET <- Li.countsF
  }
  if(the.tissue == "cerebrum" & the.sex == "female")
  {
    # testing
    the.timepoint <- CeFe.time
    counts.SUBSET <- Ce.countsF
  }
  if(the.tissue == "kidney" & the.sex == "female")
  {
    # testing
    the.timepoint <- KiFe.time
    counts.SUBSET <- Ki.countsF
  }
  if(the.tissue == "ovary" & the.sex == "female")
  {
    # testing
    the.timepoint <- ReFe.time
    counts.SUBSET <- Re.countsF
  }
  if(the.tissue == "skeletalmuscle" & the.sex == "female")
  {
    # testing
    the.timepoint <- SmFe.time
    counts.SUBSET <- Sm.countsF
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  # convert back to counts with no log
  counts.SUBSET <- 2^(   counts.SUBSET )-0.25 
  
  
  
  
  
  
  
  
  
  library(limma)
  # Normalization. Perform voom normalization:
  group <- factor(the.timepoint)
  length(group) # 8 samples
  design <- model.matrix(~group)   # works = group 0 for 3 BWH samples, the rest are group 1
  y <- voom(counts.SUBSET,design,plot=TRUE)
  # plotMDS(y,xlim=c(-2.5,2.5))
  plotMDS(y,xlim=c(-5,5))
  
  # Linear model fitting and differential expression analysis. Fit linear models to genes and assess differential expression using the eBayes moderated t statistic. 
  fit <- eBayes(lmFit(y,design), trend=TRUE)  
  # eBayes assumes 0.01 are DEG, so here 
  class(counts.SUBSET) # matrix
  length(counts.SUBSET[,1]) # is 
  # strict
  outTableF <-topTable(fit,coef=2,number=10000,p.value=0.01,lfc=1.2,adjust.method="BH") 
  outTableF <-topTable(fit,coef=2,number=10000,p.value=0.01,lfc=(log2(1.5)),adjust.method="none") 
  # View(outTableF)
  
  # 2^0.56
  # log2(1.5)
  # change so that the focus is moribund
  outTableF$logFC <- outTableF$logFC * -1
  outTableF["GHR",]
  
  
  # save output
  outTableF$gene <- row.names(outTableF) 
  # move gene column to position 1
  outTableF <- outTableF[,c(which(colnames(outTableF)=="gene"),which(colnames(outTableF)!="gene"))]
  # the.file.name <- paste("./output/",target.tissue,"_DEG_by_sex.tsv",sep="")
  the.file.name <- paste(the.sex,"_",the.tissue,"_DEG_by_sex.tsv",sep="")
  write.table(outTableF, the.file.name, row.names=F, col.names=T, sep="\t", quote=F)
  
  
  
  # GHR drops dramatically in the moribund!
  counts.SUBSET["GHR",]
  
  
} 
# END OF BIG LOOP