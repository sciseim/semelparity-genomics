rm(list=ls()) 

setwd("/Users/sciseim/Dropbox/seimlab at NNU/antechinus genome project/--data and analysis--/-molecular evolution/CNEE analysis/CNEEs 160520/")

source("drawAlign_function.R") # from PhyloAcc GitHub

### Read in tree data
treeData <- prepare_data(tree_path = "neut_ver1_1_FINAL.mod", species_name = "model1_species_names.txt", common_name = "mammalname.txt")

### Generate evolutionary pattern and sequence alignment for one element from PhyloAcc outputs 
#### read in BF scores as well as marginal log likelihood under null, accelerated and full model #### 
score <- read.table("50264.sig.block.model1/model1/model1_elem_lik.txt", header=T)
## order score by BF1
score <- score[order(-score$logBF1),]
#### read in posteriors of substitution rates and latent conservation states ####
postZ <- read.table("50264.sig.block.model1/model1/model1_rate_postZ_M2.txt", header=T, check.names = F) 

# Select an ratite-accelerated element (e.g. logBF2 > 1 and large BF1) and plot: 
sel <- which(score$logBF2 > 1)[1] # select the element with largest BF1 and BF2 > 0
lk = score[sel, ]
k = score[sel, 1] # get the No. of the selected element
targets = c("Gracilinanus_agilis","Antechinus_flavipes","Antechinus_argentus","Antechinus_arktos") # target species
Z = unlist(postZ[postZ$No. == k, -1]) # get the posteriors of conservation states
tit = paste("logBF1:", round(lk$logBF1), "logBF2:",round(lk$logBF2), "  ") # use BF scores and posterior substitution rates as title
plotZPost(Z, treeData, target_species=targets, tit=tit, offset=5,cex.score = 2) # offset= 6 indicates the posterior of Z start from 7th column

# ```plotAlign``` function will show the sequence alignment for an element as a heatmap. 
# It will need 1) a bed file and 2) sequence alignments. 
bed <- read.delim("input.bed", header=F)
fasta <- read.alignment(file = "input.fa", format = "fasta")  
align <- as.matrix(fasta)
align <- align[treeData$tree$tip.label,]  # reorder species in the alignment to be the same as tips of the tree. The name of the species in the alignment file has to the same as in the tree!


# To plot the substitutions (as well as indels and unknown base pairs 'N') of the kth element across species, 
plotAlign(k, align, bed, treeData, target_species=targets)




# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# expect these in the output
logBF1.cutoff = 5
sig.output.df <- score[(which(score$logBF1>logBF1.cutoff & score$logBF2 > 1)),] # 31
# View(sig.output.df)
library("openxlsx")
## Create a blank workbook
wb <- createWorkbook()
addWorksheet(wb, sheetName = "PhyloAcc")
writeData(wb, "PhyloAcc", sig.output.df)
## Save workbook to working directory
#@ saveWorkbook(wb, file = "./PhyloAcc_sig_CNEEs.xlsx", overwrite = TRUE) 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
length(sig.output.df$No.) # 31









### Expected number of independent accelerations and accelerated species per element
# Next, we compute the expected number of independent accelerations and accelerated species within target species. We will first select all the ratite-accelerated elements (in this simulation, all elements are ratite accelerated and expectedly, the algorithm selects almost all the simulated elements), and then compute and get the histogram of expected number of independent accelerations within clades of target species as well as accelerated target species for each element.
sel <- score$No.[which(score$logBF1>logBF1.cutoff & score$logBF2 > 1)] # select the ratites accelerated elements

# e.g. [1]  7999 20462 29679 23007  8108 29647 13979 32001  4294
score[sel,]









### Evolutionary patterns over all elements
# Finally, we can get the overall posterior probability of acceleration on each branch for all the input elements. The higher means that branch is more likely to be accelerated across the whole genome (or all inputs). 
topZ = postZ[postZ$No. %in% sel, seq(10, ncol(postZ), by=4)] # get the posterior of Z==2 (being accelerated) for selected elements
colnames(topZ) <- sapply(colnames(topZ), function(x) strsplit(x, "_")[[1]][1])
plotZPost_all(treeData, topZ, targets) # treeData is from first step; topZ are ratite-accelerated elements selected in the previous section; targets are target species (here, ratites)


dim(topZ)
pdf("overall posterior probability.pdf")
plotZPost_all(treeData, topZ, targets) # treeData is from first step; topZ are ratite-accelerated elements selected in the previous section; targets are target species (here, ratites)
dev.off()

# In this example, all input elements are accelerated only in ratite branches, and our model detected the pattern correctly as the average acceleration probability is nearly 1 for ratites (shown as red) but it is nearly 0 for non accelerated branches. But for the real data, since all patterns are mixed together, the average posterior probability of acceleration will not be concentrated at 1.


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


# Select an ratite-accelerated element (e.g. logBF2 > 1 and large BF1) and plot: 
# which(score$logBF2 > 1)[2] # test ... should be 3=1.283552300
# sel <- which(score$logBF2 > 1)[1] # select the element with largest BF1 and BF2 > 0
# e.g. sel is 2 here
# sel.all <- sel # 31 elements
# sel <- sel.all[1] # get number 1 here -- 7999 =  CNEE23102
# No.	ID	loglik_Null	loglik_Acc	loglik_Full	logBF1	logBF2	loglik_Max_M0	loglik_Max_M1	loglik_Max_M2
# 7999	CNEE23102	-831.47486	-788.70928	-789.80805	42.765583	1.098775	-823.11249	-782.83825	-782.83846


#  sel <- which(score$logBF2 > 1)[1] # select the element with largest BF1 and BF2 > 0
# e.g. 2
#  lk = score[sel, ]
# entire row
#  k = score[sel, 1] # get the No. of the selected element

# custom
# target.CNEE <- "CNEE23102" # dummy
for(i in sig.output.df$ID)
  {
  target.CNEE <- i
  cat(target.CNEE)
# //////////////////////////

target.score.info <- score[ which(score$ID == target.CNEE) , ]
lk = target.score.info
k = target.score.info$No # get the No. of the selected element

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# plot PP of particular CNEE!
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
targets = c("Gracilinanus_agilis","Antechinus_flavipes","Antechinus_argentus","Antechinus_arktos") # target species
Z = unlist(postZ[postZ$No. == k, -1]) # get the posteriors of conservation states
tit = paste("logBF1:", round(lk$logBF1), "logBF2:",round(lk$logBF2), "  ") # use BF scores and posterior substitution rates as title

the.filename <- paste(target.CNEE,"_PP_plot.pdf",sep="")
pdf(the.filename)
plotZPost(Z, treeData, target_species=targets, tit=tit, offset=5,cex.score = 2) # offset= 6 indicates the posterior of Z start from 7th column
dev.off()



# ```plotAlign``` function will show the sequence alignment for an element as a heatmap. 
# It will need 1) a bed file and 2) sequence alignments. 
# bed <- read.delim("input.bed", header=F)
# fasta <- read.alignment(file = "input.fa", format = "fasta")  
# align <- as.matrix(fasta)
# align <- align[treeData$tree$tip.label,]  # reorder species in the alignment to be the same as tips of the tree. The name of the species in the alignment file has to the same as in the tree!


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# To plot the substitutions (as well as indels and unknown base pairs 'N') of the kth element across species, 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
the.filename <- paste(target.CNEE,"_substitutions_plot.pdf",sep="")
pdf(the.filename)
plotAlign(k, align, bed, treeData, target_species=targets)
dev.off()

}
# //////////////////////////
