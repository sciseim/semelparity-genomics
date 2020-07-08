
# create output dir
# system("mkdir -p ./output")



# the.tissues <- unique(curated.dataset$tissue)
# the.tissues <- c("cerebrum","gastrointestine","kidney","liver","skeletalmuscle","spleen","stomach","sternalgland","adrenalgland")




# loop through

# "cerebrum","gastrointestine","kidney","liver","skeletalmuscle","spleen","stomach","sternalgland",
# the.tissues <- c("adrenalgland")
# target.tissue <- i



# for (i in the.tissues) {
# target.tissue <- i
#  target.tissue <- "liver"






# ######################################################################
#
#             LIMMA AND FINAL FILTERING!
#
# ######################################################################

# voom (and limma more generally) require replicates. The whole purpose of voom is to estimate the mean-variance relationship. It would work if you had any replicates at all in any of groups
# This approach is more flexible in adapting automatically to gene-specific variability in 
# RNA-Seq data than the edgeR algorithm, and has proved successful on some 
# high-variability datasets = post by Gordon Smyth
# check if GHRLOS passes the filter with two vector samples


curated.dataset.SUBSET <- curated.dataset[which(curated.dataset$tissue == target.tissue), ]
length(curated.dataset.SUBSET$samplename) # 14
sample.names.of.interest <- curated.dataset.SUBSET$samplename
counts.SUBSET <- count.table[,sample.names.of.interest]

group <- factor(curated.dataset.SUBSET$sex)
length(group) # 14 samples


count.cut.off <- 10
nosamples <- length((counts.SUBSET[1,]))
counts.SUBSET <- counts.SUBSET[rowSums(counts)>3,] # was 3 ... not that strict!
class(counts.SUBSET) # df
length(counts.SUBSET)/nosamples # 7,608 genes

library(limma)
# Normalization. Perform voom normalization:
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
outTableF <-topTable(fit,coef=2,number=10000,p.value=0.01,lfc=1.5,adjust.method="BH") 
outTableF <-topTable(fit,coef=2,number=10000,p.value=0.05,lfc=0.56,adjust.method="fdr") 
outTableF <-topTable(fit,coef=2,number=10000,p.value=0.25,lfc=2,adjust.method="fdr") 


# NOT VERY STRICT BUT USEFUL TO EXAMINE
# “holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”

outTableF <-topTable(fit,coef=2,number=10000,p.value=0.05,lfc=2,adjust.method="none") 

# outTableF <-topTable(fit,coef=2,number=10000,p.value=0.05,lfc=2,adjust.method="bonferroni") 

# outTableF <-topTable(fit,coef=2,number=10000,p.value=0.05,lfc=1.5,adjust.method="BH") 
 outTableF <-topTable(fit,coef=2,number=10000,p.value=0.01,lfc=1.5,adjust.method="none") 
 


# save output
outTableF$gene <- row.names(outTableF) 
# move gene column to position 1
outTableF <- outTableF[,c(which(colnames(outTableF)=="gene"),which(colnames(outTableF)!="gene"))]
the.file.name <- paste("./output/",target.tissue,"_DEG_by_sex.tsv",sep="")
write.table(outTableF, the.file.name, row.names=F, col.names=T, sep="\t", quote=F)




# heat map
keepgenes <- as.character(row.names(outTableF))
mergedDF <- counts.SUBSET[row.names(counts.SUBSET) %in% keepgenes, ]
mergedMATRIX <- mergedDF
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# get rid of blanks, NAs
# remember, we already have a matrix
class(mergedMATRIX) 
mergedMATRIX[mergedMATRIX==""] <- 0  # kill blanks
mergedMATRIX[is.na(mergedMATRIX)] <- 0 # and kill NAs 
# ... the above will choke heatmap
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Load necessary packages
library("RColorBrewer")
library("gplots")
library("devtools")
#Load latest version of heatmap.3 function
# source_url("https://raw.githubusercontent.com/obigriffith/biostar-tutorials/master/Heatmaps/heatmap.3.R")
source("heatmap.3.R")
#Define custom dist and hclust functions for use with heatmaps
mydist=function(c) {dist(c,method="euclidian")}
myclust=function(c) {hclust(c,method="average")}
cols <- rev(colorRampPalette(brewer.pal(10,"RdBu"))(256))
cols<- colorRampPalette(c("red", "white", "blue"))(256)
# try to transpose instead
the.file.name <- paste("./output/",target.tissue,"_DEG_heatmap.pdf",sep="")
sex <- curated.dataset.SUBSET$sexcolor
clab <- (cbind(sex))
clab <- as.matrix(clab)

pdf(the.file.name, height=15, width=10)
h <- heatmap.3((mergedMATRIX), col=cols, scale="row", trace="none",density.info="none", dendrogram="none",Colv=FALSE, key=TRUE,cexRow=0.5,breaks=seq(-3,3,6/256) , hclustfun=myclust, distfun=mydist,ColSideColors=clab,ColSideColorsSize=1)
# h <- heatmap.3((mergedMATRIX), col=cols, scale="row", trace="none",density.info="none", dendrogram="none",Colv=FALSE, key=TRUE,cexRow=0.5, hclustfun=myclust, distfun=mydist)
# breaks=seq(-3,3,6/256) ,
print(h)
dev.off()


# } # END OF TISSUE LOOP


# plot a gene of interest
geneofinterest <- "NOG"
boxplot(log(as.matrix(counts.SUBSET)[geneofinterest,])~factor(group)) # log-scale (CPM)
# boxplot((as.matrix(counts.SUBSET)[geneofinterest,])~factor(group)) # CPM 


