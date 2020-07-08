# 7B-GO and pathway enrichment analysisV4.R

setwd("~/Downloads/MIT-msigdb/")

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# load all normalised counts
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
counts.LOG <- readRDS("normalized_counts.rds")
# load("input/HPAdataFinal.Robj")
idAll <- rownames(counts.LOG)

# create output folder
system("mkdir -p ./output")
system("mkdir -p ./output/tables")
system("mkdir -p ./output/genes")

source("MIT-msigdb/FUN.hyperGTest.R")



# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# decode on gene sets to interrogates 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# gmt.databases <- c("immuno.sig","hallmark","KEGG","reactome","cis.motifs","GO.BP","GO.CC","CGP","custom")

# gmt.databases <- c("hallmark") # 50 genes. Run this first to check that all is good
# gmt.databases <- c("hallmark","KEGG","reactome","cis.motifs","GO.BP","GO.CC","CGP","custom","immuno.sig")
gmt.databases <- c("KEGG")



# gmt.databases <- c("reactome","cis.motifs","GO.BP","GO.CC","CGP","custom","immuno.sig")


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# load input file
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
files <- list.files(pattern = "\\.txt$")
# loop 

# the.database <- "hallmark"  # testing
# the.filename <- thefile
  #  the.filename <- "moribund-skeletalmuscle-male.txt"
#  the.filename <- "CORT-skeletalmuscle.txt"
 #  the.filename <- "CORT-liver.txt"
   the.filename <- "CORT-kidney.txt"
  # the.filename <- "CORT-cerebrum.txt"
 # the.filename <- "moribund-liver-male.txt"  
  #  
 
 # all.curated.gene.sets
 # START OF BIG LOOP
#@ for (thefile in files) {
   # the.database <- "hallmark"  # testing
# @   the.filename <- thefile
 
   # the.filename <- "moribund-liver-male.txt"
   # the.filename <- "moribund-skeletalmuscle-male.txt"
   # the.filename <- "moribund-kidney-male.txt"
   #  the.filename <- "moribund-cerebrum-male.txt"
   #  the.filename <- "moribund-testis-male.txt"
   
   # the.filename <- "CORT-skeletalmuscle.txt"
  #  the.filename <- "CORT-liver.txt"
   #  the.filename <- "CORT-kidney.txt"
   # the.filename <- "CORT-cerebrum.txt"

# the.filename <- "moribund-skeletalmuscle-female.txt"
# the.filename <- "moribund-kidney-female.txt"
# the.filename <- "moribund-cerebrum-female.txt"
# the.filename <- "moribund-ovary-female.txt"



   
   
   
   
   
   
# for each sample 
  DEG.df <- read.delim2(the.filename,header=FALSE) 
  DEG.df$V1 <- as.character(DEG.df$V1) 
  DEG.df$V2 <- as.character(DEG.df$V2) 
  targetgenes <- DEG.df$V1
  output.name <- gsub(".txt","",the.filename)
  output.name <- paste(output.name,"-all",sep="")
  remove(results)
  remove(enrichment)
  tryCatch(source("SUB-enrich.R"))
  # if we want to subset by direction of gene expression
  DEG.df.UP <- DEG.df[   which(DEG.df$V2 == "UP"), ]
  targetgenes <- DEG.df.UP$V1
  output.name <- gsub(".txt","",the.filename)
  output.name <- paste(output.name,"-UP",sep="")
  remove(results)
  remove(enrichment)
  tryCatch(source("SUB-enrich.R"))
  DEG.df.DOWN <- DEG.df[   which(DEG.df$V2 == "DOWN"), ]
  targetgenes <- DEG.df.DOWN$V1
  output.name <- gsub(".txt","",the.filename)
  output.name <- paste(output.name,"-DOWN",sep="")
  remove(results)
  remove(enrichment)
  tryCatch(source("SUB-enrich.R"))



  
#@  } # END OF FILE LOOP

