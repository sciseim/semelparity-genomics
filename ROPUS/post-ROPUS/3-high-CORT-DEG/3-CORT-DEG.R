rm(list=ls()) # reset R  
setwd("~/Downloads/testing/HERE/")
owd <- getwd()

library(gtools) # for SUB.CORT-absFC-v2.R

library(reshape2)
source("SUB-summarySE.R")

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# load all normalised counts
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
counts.LOG <- readRDS("normalized_counts.rds")


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




run.bootstrap <- "no" 
use.non.log.values <- "no" # use log2(TMM-normalized+0.25) 

# ############################### 
# run analysis for each tissue
# ############################### 
# combine male and female
the.tissue <- "skeletalmuscle"
Sm.df <- cbind(Sm.countsM,Sm.countsF)
names(Sm.df)
low.CORT.F <- subset(Sm.df, select=c(FN21818Sm,FN21818Sm2,FN21818Sm3,FB22918Sm,FB25918Sm,FB21018Sm2))
low.CORT.M <- subset(Sm.df, select=c(MN21818Sm,MN21818Sm))
high.CORT <- subset(Sm.df, select=-c(MN21818Sm,FN21818Sm,FN21818Sm2,FN21818Sm3,FB22918Sm,FB25918Sm,FB21018Sm2))
# source("SUB.CORT-highMinusLow.R")
source("SUB.CORT-absFC-v2.R")

target.gene <- "HLA-DPA1"
Sm.df[target.gene,]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
the.tissue <- "liver"
Li.df <- cbind(Li.countsM,Li.countsF)

target.gene <- "ZNF451"
Li.df[target.gene,]


low.CORT["PCK1",]
high.CORT["PCK1",]

2^(Li.df["KCNJ15",] -0.25 )
# log2(2)+0.25
# 2^(1.25-0.25)

names(Li.df)
low.CORT.F <- subset(Li.df, select=c(FN21818Li,FN21818Li2,FB22918Li,FB25918Li,FB021018Li2))
low.CORT.M <- subset(Li.df, select=c(MN21818Li,MN21818Li))
high.CORT <- subset(Li.df, select=-c(MN21818Li,FN21818Li,FN21818Li2,FB22918Li,FB25918Li,FB021018Li2))
# source("SUB.CORT-highMinusLow.R")
source("SUB.CORT-absFC-v2.R")
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
the.tissue <- "kidney"
Ki.df <- cbind(Ki.countsM,Ki.countsF)
names(Ki.df)
low.CORT.F <- subset(Ki.df, select=c(FN21818Ki,FN21818Ki2,FN21818Ki3,FB22918Ki,FB25918Ki,FB21018Ki2))
low.CORT.M <- subset(Ki.df, select=c(MN21818Ki,MN21818Ki))
high.CORT <- subset(Ki.df, select=-c(MN21818Ki,FN21818Ki,FN21818Ki2,FN21818Ki3,FB22918Ki,FB25918Ki,FB21018Ki2))
# source("SUB.CORT-highMinusLow.R")
source("SUB.CORT-absFC-v2.R")

target.gene <- "SLC25A30"
Ki.df[target.gene,]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
the.tissue <- "cerebrum"
Ce.df <- cbind(Ce.countsM,Ce.countsF)
names(Ce.df)
low.CORT.F <- subset(Ce.df, select=c(FN21818Ce,FN21818Ce2,FB22918Ce,FB25918Ce,FB021018Ce2))
low.CORT.M <- subset(Ce.df, select=c(MN21818Ce,MN21818Ce))
high.CORT <- subset(Ce.df, select=-c(MN21818Ce,FN21818Ce,FN21818Ce2,FB22918Ce,FB25918Ce,FB021018Ce2))
# source("SUB.CORT-highMinusLow.R")
source("SUB.CORT-absFC-v2.R")

target.gene <- "NDUFA4L2"
Ce.df[target.gene,]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


