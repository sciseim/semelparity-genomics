rm(list=ls()) # reset R  
library(reshape2)


setwd("~/Dropbox/seimlab at NNU/antechinus genome project/--data and analysis--/---late_2019---/ROPUS_data/flavipes_RNAseq/FCM/") # OSX
owd <- getwd()

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# load all normalised counts
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
counts.LOG <- readRDS("../normalized_counts.rds")


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# SET VARIABLES
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
remove(targetgene)
targetgene <- "ALB"
the.sex <- "female"



# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# KEEP ONLY MALE OR FEMALE SAMPLES
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# rationale: keep only samples samples used in the FCM!
counts.by.samling.date.female <- counts.LOG[,c("FN21818Li","FN21818Li2","FB22918Li","FB25918Li","FB021018Li2",
                                                 "FN21818Sm","FN21818Sm2","FN21818Sm3","FB22918Sm","FB25918Sm","FB21018Sm2",
                                                 "FN21818Ki","FN21818Ki2","FN21818Ki3","FB22918Ki","FB25918Ki","FB21018Ki2",
                                                 "FN21818Ce","FN21818Ce2","FB22918Ce","FB25918Ce","FB021018Ce2",
                                                 "FN21818Ov","FN21818Ov2","FB22918Ov","FB25918Ov","FB21018Ov2")]
  counts.by.samling.date.TEMP <- counts.by.samling.date.female
  # removed FN120517Ov
  names(counts.by.samling.date.TEMP)


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# SUBSET SO THAT THE *SINGLE* GENE OF INTEREST ONLY IS RETAINED
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# counts for ALL samples counts.df2
# want non-log2 counts for error calculation
# subset by genes of interest
counts.by.samling.date.TEMP.noLog <- 2^(  ( counts.by.samling.date.TEMP[which(row.names(counts.by.samling.date.TEMP) %in% targetgene),]) )-0.25
# class(counts.df2) # df
# counts.df2 <-   counts.by.samling.date.TEMP[which(row.names(counts.by.samling.date.TEMP) %in% targetgene),] 
class(counts.by.samling.date.TEMP.noLog) # df


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# SET TIMES
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# female times
LiFe.time=c("21818","21818","22918","25918","021018")
SmFe.time=c("21818","21818","21818","22918","25918","021018")
KiFe.time=c("21818","21818","21818","22918","25918","021018")
CeFe.time=c("21818","21818","22918","25918","021018")
ReFe.time=c("21818","21818","22918","25918","021018")


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# KEEP ONLY MALE OR FEMALE SAMPLES FOR EACH FIGURE ...
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# ... for FEMALES
# assign into one df per tissue
  Li.counts <- counts.by.samling.date.TEMP.noLog[,  grep("Li", names(counts.by.samling.date.TEMP.noLog) )  ]
  Sm.counts <- counts.by.samling.date.TEMP.noLog[,  grep("Sm|Sk", names(counts.by.samling.date.TEMP.noLog) )  ]
  Ki.counts <- counts.by.samling.date.TEMP.noLog[,  grep("Ki", names(counts.by.samling.date.TEMP.noLog) )  ]
  Ce.counts <- counts.by.samling.date.TEMP.noLog[,  grep("Ce", names(counts.by.samling.date.TEMP.noLog) )  ]
  Re.counts <- counts.by.samling.date.TEMP.noLog[,  grep("Ov", names(counts.by.samling.date.TEMP.noLog) )  ]
  
  
  
  # ~~~~~~~~~
  # liver
  counts.by.samling.date.subset.gene <- Li.counts # assign
  time <- LiFe.time
  remove(exprs_with_time)
  source("SUB-expression-by-time.R") # so, average expression at a time point
  remove(error.per.sample)
  source("SUB-pull-out-error.R") # pull out the error per time-point
  # assign to variable for later
  exprs_with_time.liver <- exprs_with_time
  error.per.sample.liver <- error.per.sample
  # ~~~~~~~~~
  # ~~~~~~~~~
  # skeletal muscle
  counts.by.samling.date.subset.gene <- Sm.counts # assign
  time <- SmFe.time
  remove(exprs_with_time)
  source("SUB-expression-by-time.R") # so, average expression at a time point
  remove(error.per.sample)
  source("SUB-pull-out-error.R") # pull out the error per time-point
  # assign to variable for later
  exprs_with_time.skeletalmuscle <- exprs_with_time
  error.per.sample.skeletalmuscle <- error.per.sample
  # ~~~~~~~~~
  # ~~~~~~~~~
  # kidney
  counts.by.samling.date.subset.gene <- Ki.counts # assign
  time <- KiFe.time
  remove(exprs_with_time)
  source("SUB-expression-by-time.R") # so, average expression at a time point
  remove(error.per.sample)
  source("SUB-pull-out-error.R") # pull out the error per time-point
  # assign to variable for later
  exprs_with_time.kidney <- exprs_with_time
  error.per.sample.kidney <- error.per.sample
  # ~~~~~~~~~
  # ~~~~~~~~~
  # cerebrum
  counts.by.samling.date.subset.gene <- Ce.counts # assign
  time <- CeFe.time
  remove(exprs_with_time)
  source("SUB-expression-by-time.R") # so, average expression at a time point
  remove(error.per.sample)
  source("SUB-pull-out-error.R") # pull out the error per time-point
  # assign to variable for later
  exprs_with_time.cerebrum <- exprs_with_time
  error.per.sample.cerebrum <- error.per.sample
  # ~~~~~~~~~
  # ~~~~~~~~~
  # reproductive
  counts.by.samling.date.subset.gene <- Re.counts # assign
  time <- ReFe.time
  remove(exprs_with_time)
  source("SUB-expression-by-time.R") # so, average expression at a time point
  remove(error.per.sample)
  source("SUB-pull-out-error.R") # pull out the error per time-point
  # assign to variable for later
  exprs_with_time.repro <- exprs_with_time
  error.per.sample.repro <- error.per.sample
  # ~~~~~~~~~
  
  
  
  
  
  # below works for any
  # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  # combine back into DFs
  # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  the.combined.df <- as.data.frame ( cbind(exprs_with_time.liver,exprs_with_time.skeletalmuscle,exprs_with_time.kidney,exprs_with_time.cerebrum,exprs_with_time.repro) )
  # row.names(the.combined.df) <- targetgene
  melted.df <- melt(the.combined.df)
  head(melted.df) # correct order and by time
  
  # add these once combined.df is melted 
  the.combined.errors <-  c(error.per.sample.liver,error.per.sample.skeletalmuscle,error.per.sample.kidney,error.per.sample.cerebrum,error.per.sample.repro)  # combine errors 
  length(melted.df[,1]) # 25
  length(the.combined.errors) # 25
  melted.df$error <- the.combined.errors
  
  melted.df$error[is.na(melted.df$error)] <- 0 # change any NA to 0
  melted.df$error <- as.numeric(melted.df$error)
  
  # assign tissuetype column to melted column ... remember TIMEPOINTS
  temp.df <- data.frame(tissuetype = c(
    rep("liver",length(exprs_with_time.liver)), 
    rep("skeletalmuscle",length(exprs_with_time.skeletalmuscle)),
    rep("kidney",length(exprs_with_time.kidney)),
    rep("cerebrum",length(exprs_with_time.cerebrum)),
    rep("repro",length(exprs_with_time.repro))
  ))
  melted.df$tissuetype <- temp.df$tissuetype
  
  # assign tissuetype column to melted column
  if(the.sex == "male")
  {
    temp.df <- data.frame(tissuecolor = c(
      rep("#e6aa00",length(exprs_with_time.liver)), 
      rep("#35a02e",length(exprs_with_time.skeletalmuscle)),
      rep("#646464",length(exprs_with_time.kidney)),
      rep("#a5751a",length(exprs_with_time.cerebrum)),
      rep("#a6cee3",length(exprs_with_time.repro))
    ))
    melted.df$tissuecolor <- temp.df$tissuecolor
  }
  # assign tissuetype column to melted column
  if(the.sex == "female")
  {
    temp.df <- data.frame(tissuecolor = c(
      rep("#e6aa00",length(exprs_with_time.liver)), 
      rep("#35a02e",length(exprs_with_time.skeletalmuscle)),
      rep("#646464",length(exprs_with_time.kidney)),
      rep("#a5751a",length(exprs_with_time.cerebrum)),
      rep("#ff7de3",length(exprs_with_time.repro))
    ))
    melted.df$tissuecolor <- temp.df$tissuecolor
  }
  
  
  # add times
  if(the.sex == "male")
  {melted.df$time <- c(unique(LiMa.time),unique(SmMa.time),unique(KiMa.time),unique(CeMa.time),unique(ReMa.time)  )}
  if(the.sex == "female")
  {melted.df$time <- c(unique(LiFe.time),unique(SmFe.time),unique(KiFe.time),unique(CeFe.time),unique(ReFe.time)  )}
  
  # let us rename the melted.df columns into something sane
  # names(melted.df)[names(melted.df) == 'variable'] <- 'samplename'
  names(melted.df)[names(melted.df) == 'value'] <- 'expression'
  class(melted.df)
  
  # If errorbars overlapp, use position_dodge to move them horizontally
  # pd <- position_dodge(0.1) # move them .05 to the left and right
  pd <- position_dodge(0.05) # move them .05 to the left and right
  
  # if we want to use log values again
  melted.df$log2 <- log2(melted.df$expression+0.25)
  melted.df$log2error <- log2(melted.df$error)
  # convert -Inf to 0 -- i.e. no error due to single replicate or perfect (right!) measurements
  melted.df$log2error[which(!is.finite(melted.df$log2error))] <- 0
  
  # LiMa.time
  # [1] "21818"  "22918"  "22918"  "25918"  "25918"  "021018" "021018" "041018"
  
  bac <- melted.df
  melted.df$error[is.na(melted.df$error)] <- 0
  
  
  # not sure if this works yet ...
  # order df by time
  # melted.df <- melted.df[order(melted.df$time),]
  # melted.df$expression <- factor(melted.df$expression , levels = melted.df$expression [order(melted.df$time)])
  melted.df$time
  names(melted.df)
  
  head(melted.df[1:6,1:6])
  
  # make R respect the time order we want
  # different for female
  
  if(the.sex == "male")
  {
    melted.df$time <- factor(melted.df$time,levels = c("21818" , "22918" , "25918" , "021018" ,"041018"))
  }
  if(the.sex == "female")
  {
    melted.df$time <- factor(melted.df$time,levels = c("21818" , "22918" , "25918" , "021018"))
  }
  
  
# melted.df$tissuetype <- factor(melted.df$time,levels = as.character(unique(melted.df$tissuecolor)  )) 
  
  melted.df$tissuetype
  melted.df$tissue
  names(melted.df)
  
# note: using aes here so that any reordering in honoured in the colour scheme!
# also geom_line must go first so the plot points are in front
the.gene.plot <- ggplot(melted.df, aes(x=time, y=expression, group=tissuetype))  + theme_classic() +  labs(x = "", y = "Expression (CPM)") + geom_line(aes(colour = tissuetype),size = 2.5) + theme(axis.title.x = element_text(color = "black", size = 14, face = "bold"),axis.title.y = element_text(color = "black", size = 14, face = "bold"),axis.text.y = element_text(color = "black", size = 14) ,axis.text.x = element_text(color = "black", size = 14) ) + geom_errorbar(aes(ymin=expression-error, ymax=expression+error), colour="black", width=.01, position=pd) + geom_point(size = 0.1, stroke = 8, shape = 16,color="black")
  
  
  
  the.gene.plot.file.name <- paste(the.sex,"_",targetgene,".pdf",sep="")
  pdf(the.gene.plot.file.name)
  the.gene.plot
  dev.off()
  

  

  
  
# # LOG2
# melted.df
# the.gene.plot <- ggplot((melted.df), aes(x=time, y=log2, group=tissuetype))  + theme_classic() +  labs(x = "", y = "Expression log2(CPM)") + geom_line(aes(colour = tissuetype),size = 2.5) + theme(axis.title.x = element_text(color = "black", size = 14, face = "bold"),axis.title.y = element_text(color = "black", size = 14, face = "bold"),axis.text.y = element_text(color = "black", size = 14) ,axis.text.x = element_text(color = "black", size = 14) ) + geom_errorbar(aes(ymin=log2-log2error, ymax=log2+log2error), colour="black", width=.01, position=pd) + geom_point(size = 0.1, stroke = 8, shape = 16,color="black")
# the.gene.plot.file.name <- paste(the.sex,"_",targetgene,"-LOG2.pdf",sep="")
#   pdf(the.gene.plot.file.name)
#   the.gene.plot
#   dev.off()
  
  













