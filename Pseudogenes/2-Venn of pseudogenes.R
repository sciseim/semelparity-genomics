


flavipes.pseudogenes <- read.delim2("A.flav.daiqu-pseudopipe.output.class.shifts_et_stops.idcutoff.txt",header=TRUE)
agilis.pseudogenes <- read.delim2("G.agilis.fushu-pseudopipe.output.class.shifts_et_stops.idcutoff.txt",header=TRUE)


length(flavipes.pseudogenes$gene) # 483 
length(agilis.pseudogenes$gene) # 485

length(unique(flavipes.pseudogenes$gene)) # 339
length(unique(agilis.pseudogenes$gene)) # 379

# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
# Venn of A.flavipes and G. agilis
options(stringsAsFactors=F)

# BiocManager::install("RBGL")
# BiocManager::install("graph")
# install.packages("devtools")
library(devtools);
# install_github("js229/Vennerable")
library(Vennerable)
# install.packages("colorfulVennPlot")
library(colorfulVennPlot)
require(RColorBrewer)  
#establim els colors per als plots
cols <- brewer.pal(8,"Pastel2") 

head(flavipes.pseudogenes)

class(flavipes.pseudogenes$gene) # character

my.pseudogenes <- list("flavipes" = unique((flavipes.pseudogenes$gene)),
                       "agilis" = unique((agilis.pseudogenes$gene)) )
str(my.pseudogenes)

intersect( my.pseudogenes[[1]] ,my.pseudogenes[[2]]   ) # cerebrum and liver

overlapping.genes <- intersect( my.pseudogenes[[1]] ,my.pseudogenes[[2]]   ) # cerebrum and liver
overlapping.genes.df <- flavipes.pseudogenes[ which(flavipes.pseudogenes$gene %in% overlapping.genes) ,]
names(overlapping.genes.df)
class(overlapping.genes.df)
overlapping.genes.df <- overlapping.genes.df[,c("gene","description")]

length(unique(overlapping.genes.df$gene)) # 51
overlapping.genes.df <- (unique(overlapping.genes.df)) # 51
length((overlapping.genes.df$gene)) # 51

write.table(overlapping.genes.df, "overlapping.genes.tsv", sep="\t",row.names = FALSE,append=FALSE,col.names = TRUE,quote = FALSE)


# type


# only works for 4 
my.pseudogenes <- Venn(my.pseudogenes)
pdf("Venn_Aflavipes_Gagilis-pseudogenes.pdf")
plot(my.pseudogenes, type= "circles",doWeights=F)
dev.off()
# VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV

?Vennerable