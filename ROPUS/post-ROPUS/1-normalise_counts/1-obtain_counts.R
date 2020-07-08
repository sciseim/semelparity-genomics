# reset R  
rm(list=ls())





################################################################################################################################################################################################################


setwd("~/Dropbox/seimlab at NNU/antechinus genome project/--data and analysis--/---late_2019---/ROPUS_data/flavipes_RNAseq/CDS-count/") # OSX
# setwd("~/Downloads/HOMEWORK101219/CDS-count/") # OSX


#
# getwd()

# ######################################################################
#
#             LOAD sample
#
# ######################################################################

samplenames <- list.files(pattern = ".")
# list.files
# 145 files. One excluded since no reads!

class(samplenames)
length(samplenames) # 123
# get rid of summary
# samplenames <- samplenames[!samplenames %in% ".summary"]






# ######################################################################
#
#             REMOVE BAD SAMPLES
#
# ######################################################################
# use multiQC (will parse the .summary files)
# multiqc . 

# remove.these <- c(
#   "Antechinus_spleen_male_B_MB041018Sp",
#   "Antechinus_gastrointestinal_male_B_MB22918GIT2",
#   "Antechinus_spleen_female_B_FB28918Sp",
#   "Antechinus_spleen_male_B_MB22918Sp",
#   "Antechinus_stomach_female_NB_FN21818St3",
#   "Antechinus_gastrointestinal_male_NB_MN21818GIT",
#   "Antechinus_prostate_male_B_MB021018Pr",
#   "Antechinus_stomach_male_B_MB25918St",
#   "Antechinus_kidney_female_NB_FN21818Ki3",
#   "Antechinus_cerebrum_male_B_MB25918Ce",
#   "Antechinus_kidney_female_B_FB22918Ki",
#   "Antechinus_kidney_male_B_MB25918Ki",
#   "Antechinus_prostate_male_B_MB22918Pr",
#   "Antechinus_kidney_female_B_FB21018Ki2",
#   "Antechinus_kidney_male_B_MB021018Ki2",
#   "Antechinus_spleen_male_NB_MN290817Sp2",
#   "Antechinus_kidney_male_B_MB041018Ki",
#   "Antechinus_kidney_male_NB_MN21818Ki",
#   "Antechinus_liver_male_NB_MN21818Li",
#   "Antechinus_kidney_male_B_MB25918Ki2",
#   "Antechinus_liver_female_NB_FN21818Li3",
#   "Antechinus_kidney_female_B_FB28918Ki",
#   "Antechinus_liver_female_B_FB28918Li",
#   "Antechinus_kidney_male_B_MB22918Ki2",
#   "Antechinus_skeletalmuscle_female_B_FB28918Sm",
#   "Antechinus_skeletalmuscle_male_B_MB021018Sm",
#   "Antechinus_prostate_male_B_MB021018Pr2")
# samplenames <- subset(samplenames, !(samplenames %in% remove.these))
# 
# class(samplenames)
# length(samplenames) # 118
# 
# better to also visualise... (see box plots below)


remove.these <- c("Antechinus_gastrointestinal_male_B_MB22918GIT2",
                  "Antechinus_prostate_male_B_MB021018Pr",
                  "Antechinus_prostate_male_B_MB021018Pr2",
                  "Antechinus_prostate_male_B_MB021018Pr2",
                  "Antechinus_spleen_female_B_FB28918Sp",
                  "Antechinus_spleen_male_B_MB041018Sp",
                  "Antechinus_spleen_male_B_MB22918Sp",
                  "Antechinus_stomach_female_B_FB220918St",
                  "Antechinus_stomach_female_NB_FN21818St3",
                  "Antechinus_stomach_male_B_MB25918St",
                  "Antechinus_gastrointestinal_male_NB_MN21818GIT",
                  "Antechinus_prostate_male_B_MB25918Pr",
                  "Antechinus_adrenalgland_female_NB_FN21818Ad3",
                  "Antechinus_cerebrum_female_NB_FN21818Ce3",
                  "Antechinus_liver_female_NB_FN21818Li3",
                  "Antechinus_ovary_female_NB_FN21818Ov3",
                  "Antechinus_spleen_female_NB_FN21818Sp3",
                  "Antechinus_cerebrum_male_B_MB25918Ce"
)




# added 201219
# Antechinus_adrenalgland_female_NB_FN21818Ad3
# Antechinus_cerebrum_female_NB_FN21818Ce3
# Antechinus_liver_female_NB_FN21818Li3
# Antechinus_ovary_female_NB_FN21818Ov3
# Antechinus_spleen_female_NB_FN21818Sp3
# Antechinus_cerebrum_male_B_MB25918Ce







samplenames <- subset(samplenames, !(samplenames %in% remove.these))
# 
class(samplenames)
length(samplenames) # 107





matching <- as.data.frame(samplenames) # to make it legacy compatible with previous code
length(matching$samplenames) # 








result = list()
for (i in 1:nrow(matching)) {
  temp <- read.delim(paste("./", matching[i,1], "", sep=""), comment.char = "#", header=T)
  temp$Geneid <- gsub(".*_","",temp$Geneid)
  temp$Chr <- gsub(".*_","",temp$Chr)
  rownames(temp) <- temp[,1]
  
  head(temp)
  tail(temp)
  rownames(temp)
  
  temp <- temp[,c(1,7)]
  
  #  temp = temp[,c(1,2,6,7)] to re-order
  #  colnames(temp) = c("Ortholog.ID", "Species.ID", "Length", "Count")
  
  result[[as.character(matching[i,1])]] <- temp
}




# which gene names are common?
gene.common <- table(unlist(lapply(result, rownames)))
nosamples <- length(result)
gene.common <- names(gene.common)[gene.common==nosamples]  # change to # of samples!  
head(gene.common)
gene.common <- gsub(".mrna1","",gene.common) # # remove .mrna1 suffix




# only keep common genes in the list results (really a dfList)
result.common <- lapply(result, function(x) x[gene.common,])
head(result.common)

names(result.common)

# create the final table with counts 
count.table <- do.call(cbind.data.frame, lapply(result.common, function(x) x[,2])) # ,2 is count
rownames(count.table) <- gene.common

# pre-normalization!
pdf(file = "../pre-normalised_counts.pdf",   # The directory you want to save the file in
    width = 100, # The width of the plot in inches
    height = 20) # The height of the plot in inches
boxplot(log(count.table), las = 2)
dev.off()
# clearly some samples need to go

# better! still a couple of outlier, but not bad




head(count.table)
tail(count.table)

# kill column 1. It is stuffing up things
# count.table <- count.table[-c(1),]
head(count.table[1:3,1:4])


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# if there is a row with all NA (happens with exons)
# (count.table[1:2,])
count.table <- count.table[complete.cases(count.table), ]
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



























# ######################################################################
#
#             REMOVE OUTLIERS
#
# ######################################################################
# those ortholog sets with too high counts
# (i.e. read counts contributing to >5% of the total counts; three orthologs were removed this way) or
# too low counts (i.e. <10 counts in four or more samples) were discarded. The library sizes were the



# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# FIX SAMPLE NAMES HERE 
#
setwd("../") # go one folder down
library(openxlsx)
# load sample info
curated.dataset <- read.xlsx("sampleinfo.xlsx")
head(curated.dataset)
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@




# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# FIX NAMES, MAKE SURE SAMPLEINFO IS OK
# names(counts.LOG) <- gsub("Antechinus_","",names(counts.LOG))
# want to remove every before last underscore _
names(count.table) <- sub('.*\\_', '', names(count.table))
# only keep if a match in sampleinfo
keep.these.columns <- which (names(count.table) %in% curated.dataset$samplename)
length(keep.these.columns) # 112
count.table <- count.table[ , keep.these.columns]
length(count.table[1,]) # 112
# and vice versa!
keep.these.columns <- which (curated.dataset$samplename %in% names(count.table))
length(keep.these.columns) # 127
curated.dataset <- curated.dataset[  keep.these.columns ,]
length(curated.dataset$samplename) # 127

class(curated.dataset$sample.dated) # character
curated.dataset$sample.dated <- as.Date(curated.dataset$sample.dated, "%d/%m/%Y")
# Sort by vector name [z] then [x]
# curated.dataset[with(curated.dataset, order(tissue, as.Date(curated.dataset$sample.dated, format = "%d/%m/%Y"))), ]
curated.dataset <- curated.dataset[
  with(curated.dataset, order(as.Date(curated.dataset$sample.dated, format = "%d/%m/%Y")) ,testis ), ]
# then sort our counts in this order so that they are by date
sample.order <- curated.dataset$samplename
count.table <- count.table[,sample.order] # now sorted correctly for later
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# save output
saveRDS(count.table, "count.table_no_names_yet.rds")
count.table <- readRDS("count.table_no_names_yet.rds")

head(count.table[1:4,1:5])




# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Gene Annotation
library(biomaRt)
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# Asian biomaRt down again, so force 
mart<-useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
              dataset = "hsapiens_gene_ensembl", 
              host = "www.ensembl.org")

#building a query, requires filters, attributes and values
#listFilters shows all filters
filters <- listFilters(mart)
head(filters)
#
Listattributes <- listAttributes(mart)
head(Listattributes)
Listattributes[1:30,] # want  external_gene_name
ensembl_gene_id <- rownames(count.table) # human gene names (ensembl)
head(ensembl_gene_id)
# query time!
martoutput <- getBM(
  filters= "ensembl_gene_id", 
  attributes= c("ensembl_gene_id","external_gene_name", "description"),
  values= ensembl_gene_id,
  mart= mart)

length(martoutput$ensembl_gene_id) # 15,248
head(martoutput)

# save output
saveRDS(martoutput, "martoutput.rds")
# martoutput <- readRDS("martoutput.rds")

# Use plyr to rename column and to mergge (plyr join has more function)
# http://www.inside-r.org/packages/cran/plyr/docs/join
library(plyr) #  from CRAN
martoutput2 <- rename(martoutput, c("ensembl_gene_id" = "gene")) # rename gene, so that it matches regressiondata
head(martoutput2)

head(count.table)
outTableF2 <-  count.table # add rowname as column
outTableF2["gene"] <- rownames(count.table)


# save ...
saveRDS(count.table, "count.table_no_names_yet.rds")
count.table <- readRDS("count.table_no_names_yet.rds")

# # Now we want to merge the biomaRt information with the linear regression data!
# martoutput3 <- martoutput2
# martoutput3[, "gene"] <- rownames(outTableF)
# head(martoutput3)

# can now join!
outTableF3 <- join(outTableF2, martoutput2, type = "left", match = "first")
head(outTableF3)

# remove columns
# outTableF3[,1] <- NULL  # kill hsapiens_homolog_ensembl_gene
tail(outTableF3)
colnames(outTableF3)
# "external_gene_name" "description"        "chromosome_name"    "transcript_biotype"

head(outTableF3)
outTableF3$external_gene_name[duplicated(outTableF3$external_gene_name)] # HSPA14
outTableF3 <- outTableF3[ which(outTableF3$external_gene_name != "HSPA14") ,]
row.names(outTableF3) <- outTableF3$external_gene_name

names(outTableF3)
count.table <- subset(outTableF3, select=-c(description,external_gene_name))
head(count.table[1:2,])
count.table$gene <- row.names(outTableF3)
names(count.table)
count.table$gene

# add new gene name as rows andd drop gene columns for now
row.names(count.table) <- count.table$gene
count.table <- subset(count.table, select=-c(gene))
count.table[1:2,]
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
































# ######################################################################
#
#             NORMALISE AND FILTER
#
# ######################################################################

# Normalization time!
library("edgeR")
tmm <- calcNormFactors(count.table, method = "TMM")
lib.size <- colSums(count.table)/10^6
eff.lib.size <- lib.size * tmm
count.table.scale <- t(t(count.table) / eff.lib.size)

#
min.count = apply(count.table.scale, 1, function(x) sum(x<3))
table(min.count)
head(count.table.scale)

# get rid of rows with only zeros...
count.table.scale.no0 <- count.table.scale[ rowSums(count.table.scale)!=0, ] 
head(count.table.scale.no0)
count.table.scale <- count.table.scale.no0 
head(count.table.scale)




# replace 0s with NA
count.table.scale.nona <- count.table.scale
#@ count.table.scale.nona[count.table.scale.nona==0] = NA
count.table.scale.nona[count.table.scale.nona==0] = 0.0000001  # Better! 
# expect many to be 0 prior to hypoxia insult!
count.table.scale <- count.table.scale.nona
head(count.table.scale)


# number of 1:1 orthologs
length(count.table.scale)/nosamples # 15,263 genes LIVER


# post-normalization!
pdf(file = "normalised_counts.pdf",   # The directory you want to save the file in
    width = 100, # The width of the plot in inches
    height = 20) # The height of the plot in inches
boxplot(log(count.table.scale))
dev.off()


# Filtering. Only keep in the analysis those genes which had >10 reads per million mapped reads in at least two libraries.
counts <- count.table.scale
head(counts)























# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# NORMALISE ALL TOGETHER, OTHERWISE WGCNA ETC LATER NOT POSSIBLE... 
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#  Add small offset to each observation to avoid taking log of zero
counts.LOG <- log2(counts + 0.25)
head(counts.LOG[1:2,1:4])
# cat(names(counts.LOG),file="counts.LOG.names.txt",sep="\n")

curated.dataset$samplename
class(counts.LOG)
counts.LOG <- as.data.frame(counts.LOG)


head(count.table[1:2,1:3])
head(counts.LOG[1:2,1:3])
count.table.with.gene.column <- count.table
count.table.with.gene.column$gene <- row.names(count.table.with.gene.column)
counts.LOG.with.gene.column <- as.data.frame ( counts.LOG )
counts.LOG.with.gene.column$gene <- row.names(counts.LOG.with.gene.column)
# move gene column to position 1
count.table.with.gene.column <- count.table.with.gene.column[,c(which(colnames(count.table.with.gene.column)=="gene"),which(colnames(count.table.with.gene.column)!="gene"))]
counts.LOG.with.gene.column <- counts.LOG.with.gene.column[,c(which(colnames(counts.LOG.with.gene.column)=="gene"),which(colnames(counts.LOG.with.gene.column)!="gene"))]



# save output (for use in WGCNA, etc)
write.table(count.table.with.gene.column, "unnormalized_counts_Aflavipes.tsv", sep="\t",row.names = FALSE)
write.table(counts.LOG.with.gene.column, "log2_counts_Aflavipes.tsv", sep="\t",row.names = FALSE)
# drop gene column before proceeding with analysis
# class(counts.LOG)
# counts.LOG <- as.data.frame(counts.LOG)
# counts.LOG <- subset(counts.LOG, select=-c(gene))
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
dim(count.table.with.gene.column)
# [1] 15246   113
dim(counts.LOG.with.gene.column)
# 15246   113


# save counts for later 
# fix later ... should save with a gene column (so stuff above)
saveRDS(count.table, "unnormalized_counts.rds")
count.table <- readRDS("unnormalized_counts.rds")
saveRDS(counts.LOG, "normalized_counts.rds")
counts.LOG <- readRDS("normalized_counts.rds")
# save sampleinfo for later 
saveRDS(curated.dataset, "curated_dataset.rds")
curated.dataset <- readRDS("curated_dataset.rds")





















# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# PCA
pcs = prcomp(t(counts.LOG), center = TRUE)
# works with from above, but fails here ... 
class(counts.LOG) # df

percentVar = round(((pcs$sdev) ^ 2 / sum((pcs$sdev) ^ 2)* 100), 2) 
## PCA Plot
library(emuR)

# ?emuR



colorinfo <- curated.dataset$tissuecolour
length(unique(curated.dataset$tissuecolour)) # 12


length(colorinfo) # 106 colours
length(pcs$x[1,]) # 106


# colorinfo <- curated.dataset$sexcolor
library("ggsci")
library("ggplot2")
library("gridExtra")
p1 <- ggplot(as.data.frame(pcs$x), aes(PC1,PC2), environment = environment()) + 
  #  xlab(makeLab(percentVar[1],1)) + ylab(makeLab(percentVar[2],2)) + ggtitle(title) + 
  # xlab(label)
  # geom_point(size = 4, aes(colour = data.info.df$Sex)) +  
  # geom_point(size = 4, aes(colour = data.info.df$breeder)) +  
  
  # colour them later!
 geom_point(size = 4, aes(colour = colorinfo)) +  
#  geom_point(size = 4) +  
  
  
  #   geom_point(size = 4, aes(colour = colorinfo)) +  
  theme(legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_text(size = 12, colour = "black", face = "bold"),
        plot.title = element_text(size = 0, face ="bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold", color = "black"),
        axis.text.y = element_text(size = 12, face = "bold", color = "black"),
        # plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")
  )
# p1 + scale_color_npg() + theme_bw() +geom_text(aes(label=curated.dataset$samplename),hjust=1, vjust=0)
# p1 + scale_color_manual(values = colorinfo) + theme_bw() +geom_text(aes(label=curated.dataset$samplename),hjust=1, vjust=0)








library(RColorBrewer)
library(ggthemes)
# Define the number of colors you want
nb.cols <- length(unique(curated.dataset$tissue))
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)



# save as large PDF
pdf(file = "PCA.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches
# p1 + scale_fill_manual(values = mycolors) + theme_bw() +geom_text(aes(label=curated.dataset$samplename),hjust=1, vjust=0) +  theme_tufte() + geom_rangeframe() 

# p1 + scale_color_manual(values=(colorinfo )) + theme_bw() +geom_text(aes(label=curated.dataset$samplename),hjust=1, vjust=0) +  theme_tufte() + geom_rangeframe() 

mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))
p1 + scale_color_manual(values = mycolors) + theme_bw() +geom_text(aes(label=curated.dataset$samplename),hjust=1, vjust=0) +  theme_tufte() + geom_rangeframe() 
dev.off()

# save as large PDF
pdf(file = "PCA-nolabels.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches
# p1 + scale_fill_manual(values = mycolors) + theme_bw() +geom_text(aes(label=curated.dataset$samplename),hjust=1, vjust=0) +  theme_tufte() + geom_rangeframe() 
# p1 + scale_color_manual(values=(colorinfo )) + theme_bw() +geom_text(aes(label=curated.dataset$samplename),hjust=1, vjust=0) +  theme_tufte() + geom_rangeframe() 
mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 6))
p1 + scale_color_manual(values = mycolors) + theme_bw()  +  theme_tufte() + geom_rangeframe() 
dev.off()



# misc testing
g2 <- ggplot(as.data.frame(pcs$x), aes(PC1,PC2)) +
  geom_point(aes(color = colorinfo))
g2 + scale_color_manual(values = colorinfo) + theme_bw() +geom_text(aes(label=curated.dataset$samplename),hjust=1, vjust=0) +  theme_tufte() + geom_rangeframe() 

g2 + scale_color_manual(values = as.factor(colorinfo)) + theme_bw() +geom_text(aes(label=curated.dataset$samplename),hjust=1, vjust=0) +  theme_tufte() + geom_rangeframe() 







# save as large PDF
pdf(file = "PCA-by-sex.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches
# p1 + scale_fill_manual(values = mycolors) + theme_bw() +geom_text(aes(label=curated.dataset$samplename),hjust=1, vjust=0) +  theme_tufte() + geom_rangeframe() 
# p1 + scale_color_manual(values=(colorinfo )) + theme_bw() +geom_text(aes(label=curated.dataset$samplename),hjust=1, vjust=0) +  theme_tufte() + geom_rangeframe() 
# by sex
colorinfo <- curated.dataset$sexcolor
colorinfo <- gsub("#f032e6","blue",colorinfo)
colorinfo <- gsub("#4363d8","pink",colorinfo)
unique(curated.dataset$sexcolor)
class(curated.dataset$sexcolor)
p1 + scale_color_manual(values= colorinfo ) + theme_bw() +geom_text(aes(label=curated.dataset$samplename),hjust=1, vjust=0) +  theme_tufte() + geom_rangeframe() 
dev.off()

# no labels
pdf(file = "PCA-by-sex-nolabels.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches
# p1 + scale_fill_manual(values = mycolors) + theme_bw() +geom_text(aes(label=curated.dataset$samplename),hjust=1, vjust=0) +  theme_tufte() + geom_rangeframe() 
# p1 + scale_color_manual(values=(colorinfo )) + theme_bw() +geom_text(aes(label=curated.dataset$samplename),hjust=1, vjust=0) +  theme_tufte() + geom_rangeframe() 
# by sex
colorinfo <- curated.dataset$sexcolor
colorinfo <- gsub("#f032e6","blue",colorinfo)
colorinfo <- gsub("#4363d8","pink",colorinfo)
unique(curated.dataset$sexcolor)
class(curated.dataset$sexcolor)
p1 + scale_color_manual(values= colorinfo ) + theme_bw()  +  theme_tufte() + geom_rangeframe() 
dev.off()

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@










# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# hierarchial clustering
class(counts.LOG)
counts.LOG <- as.matrix(counts.LOG)
select = order(rowMeans(counts.LOG), decreasing=TRUE)[1:100]
highexprgenes_counts <- counts.LOG[select,]
# heatmap with condition group as labels
colnames(highexprgenes_counts)<- curated.dataset$tissue # comment out if you want sample names
pdf(file = "High_expr_genes.heatmap.pdf",   # The directory you want to save the file in
    width = 50, # The width of the plot in inches
    height = 10) # The height of the plot in inches
heatmap(highexprgenes_counts, col=topo.colors(50), margin=c(10,6))
dev.off()

