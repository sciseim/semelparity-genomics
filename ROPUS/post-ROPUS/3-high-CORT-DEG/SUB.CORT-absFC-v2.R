# ############
# SUB
# ############


# ################################################################################################################################### 
# Method adapted from Simola and colleagues. Epigenetic (re)programming of caste-specific behavior in the ant Camponotus floridanus 
# Science 351, aac6633 (2016)
# ################################################################################################################################### 


# disable scientific notation
options(scipen=999)


# output.data <- cbind(the.gene,low.CORT.mean.F,low.CORT.mean.M,high.CORT.mean,the.fold.change.F,the.fold.change.M,the.se,the.Mann.Whitney.Utest.F,the.Mann.Whitney.Utest.M)


library(stringr)
new.df <- data.frame(gene=(character()),
                     low.CORT.F.avg=numeric(), 
                     low.CORT.M.avg=numeric(),
                     high.CORT.avg=numeric(),
                     AbsFc.high.CORT.vs.low.CORT.F=numeric(), 
                     AbsFc.high.CORT.vs.low.CORT.M=numeric(),
                     SEM=numeric(), 
                     P.value.high.CORT.vs.low.CORT=(numeric())
                   
                     ) 
                   #         P.value.high.CORT.vs.low.CORT.F=(numeric()),
                   # P.value.high.CORT.vs.low.CORT.M=(numeric())


# then for each row ...
for (i in 1:nrow(high.CORT)) 
{
  
  # testing. This gene should pass
  # which(row.names(high.CORT) == "HSD17B2")
  # i <- 1502
  
  # should fail
  # which(row.names(high.CORT) == "UBR7")
  # i <- 288
  

  
  the.row.high.CORT <- high.CORT[i,]
  the.row.low.CORT.F <- low.CORT.F[i,]
  the.row.low.CORT.M <- low.CORT.M[i,]
  
  
  if(use.non.log.values == "yes")
{
  # obtain normalised counts
  the.row.high.CORT <- 2^(  the.row.high.CORT -0.25 )
  the.row.low.CORT.F <- 2^(  the.row.low.CORT.F -0.25 )
  the.row.low.CORT.M <- 2^(  the.row.low.CORT.M -0.25 )
  }  
  
  
  # SE
  the.row.CORT.combined <- cbind(the.row.high.CORT,the.row.low.CORT.F,the.row.low.CORT.M)
  temp <- melt(the.row.CORT.combined)
  temp$group <- "combined"
  # head(temp)
  the.stats <- summarySE(temp, measurevar="value", groupvars=c("group"))
  the.se <- the.stats$se 
  the.se.cutoff <- 4.5 * the.se
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # high-CORT
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # length(the.row.high.CORT)
  temp <- melt(the.row.high.CORT)
  temp$group <- "high"
  the.stats <- summarySE(temp, measurevar="value", groupvars=c("group"))
  # high.CORT.se <- the.stats$se 
  high.CORT.mean <- mean(the.stats$value)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # low-CORT.F
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # length(the.row.low.CORT)
  temp <- melt(the.row.low.CORT.F)
  temp$group <- "low"
  the.stats <- summarySE(temp, measurevar="value", groupvars=c("group"))
  # low.CORT.se <- the.stats$se 
  low.CORT.mean.F <- mean(the.stats$value)
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # low-CORT.M
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # length(the.row.low.CORT)
  temp <- melt(the.row.low.CORT.M)
  temp$group <- "low"
  the.stats <- summarySE(temp, measurevar="value", groupvars=c("group"))
  # low.CORT.se <- the.stats$se 
  low.CORT.mean.M <- mean(the.stats$value)
  
  
  
  # high-CORT vs .low-CORT.FEMALE
  the.fold.change.F <- foldchange(  high.CORT.mean, low.CORT.mean.F   ) # so that FC are negative if down in high  CORT
  # and Mann-Whitney U-test
  the.Mann.Whitney.Utest.F <- 1.0 # assign dummy 
  the.Mann.Whitney.Utest.F <- wilcox.test(as.numeric(the.row.low.CORT.F), as.numeric(the.row.high.CORT))
  the.Mann.Whitney.Utest.F <- the.Mann.Whitney.Utest.F$p.value 
  # high-CORT vs .low-CORT.MALE
  the.fold.change.M <- foldchange(  high.CORT.mean, low.CORT.mean.M   ) # so that FC are negative if down in high  CORT
  # and Mann-Whitney U-test
  the.Mann.Whitney.Utest.M <- 1.0 # assign dummy 
  
  the.row.low.CORT.M # want to avoid ties ... will influence the Mann-Whitney
  #   the.Mann.Whitney.Utest.M <- wilcox.test(as.numeric(the.row.low.CORT.M), as.numeric(the.row.high.CORT))
  the.Mann.Whitney.Utest.M <- wilcox.test(as.numeric(the.row.low.CORT.M[,1]), as.numeric(the.row.high.CORT))
  #
  the.Mann.Whitney.Utest.M <- the.Mann.Whitney.Utest.M$p.value 
  
  
  # Science MS = combined Mann-Whitney of treatment vs control groups!
  the.Mann.Whitney.Utest.combo <- 1.0 # assign dummy 
  the.row.low.CORT.combo <- cbind(the.row.low.CORT.F,the.row.low.CORT.M)
  the.Mann.Whitney.Utest.combo <- wilcox.test(as.numeric(the.row.low.CORT.combo), as.numeric(the.row.high.CORT))
  the.Mann.Whitney.Utest.combo <- the.Mann.Whitney.Utest.combo$p.value 
  
  
  
  
  
  
  
  
  the.Mann.Whitney.Cutoff <- 0.05
  
  
  if(use.non.log.values == "yes")
  {
  normalised.count.cutoff <- 2 # NOT log2+0.25
  }
  
  normalised.count.cutoff <- log2(2+0.25)
  
  
  # and check if both fold-changes negative or both positive 
  # sign(the.fold.change.F)==sign(the.fold.change.M) == TRUE
  the.FC.test <- sign(the.fold.change.F)==sign(the.fold.change.M)
  class(the.FC.test) # logical
  # if(the.FC.test == TRUE){cat("PASS")}



  
  # ASSESS 
  if( 
    the.FC.test == TRUE &
    abs(the.fold.change.F) >= the.se.cutoff & 
    abs(the.fold.change.M) >= the.se.cutoff & 
    # the.Mann.Whitney.Utest.F <= the.Mann.Whitney.Cutoff & 
    # the.Mann.Whitney.Utest.M <= the.Mann.Whitney.Cutoff & 
    the.Mann.Whitney.Utest.combo <= the.Mann.Whitney.Cutoff & 
    abs(the.fold.change.F) != 0 & 
    abs(the.fold.change.M) != 0 &
    the.se.cutoff != 0 &
    low.CORT.mean.M >= normalised.count.cutoff &
    low.CORT.mean.F >= normalised.count.cutoff &
    high.CORT.mean >= normalised.count.cutoff
    )
  {
    # cat("PASS")
    # want to add these
    the.gene <- row.names(the.row.high.CORT)
    # low.CORT.mean
    # high.CORT.mean
    # the.fold.change
    # the.se
    # the.Mann.Whitney.Utest
 #   output.data <- cbind(the.gene,low.CORT.mean.F,low.CORT.mean.M,high.CORT.mean,the.fold.change.F,the.fold.change.M,the.se,the.Mann.Whitney.Utest.F,the.Mann.Whitney.Utest.M)
    output.data <- cbind(the.gene,low.CORT.mean.F,low.CORT.mean.M,high.CORT.mean,the.fold.change.F,the.fold.change.M,the.se,the.Mann.Whitney.Utest.combo)
    
      # add to df
    new.df <- rbind(new.df,output.data)
  } # END OF ASSESSMENT
  # 
} # END OF PER GENE LOOP



# additional fold-change cut-off limit for BOTH groups?
if(use.non.log.values == "no")
{
post.FC.cutoff <- 1.50
new.df.FC.subset <- new.df[ which( 
  abs(as.numeric(as.character(new.df$the.fold.change.F))) >= post.FC.cutoff 
  &
    abs(as.numeric(as.character(new.df$the.fold.change.M))) >= post.FC.cutoff 
), ] 
length(new.df.FC.subset[,1]) # 43
new.df <- new.df.FC.subset
}


the.filename <-paste("CORT","_",the.tissue,"-AbsFC.tsv",sep="")
write.table(new.df, the.filename, sep="\t",row.names = FALSE,append=FALSE)



















# ################################################################################################
# draw heat map
# ################################################################################################
class(the.Mann.Whitney.Utest.combo)
class(the.Mann.Whitney.Cutoff)
new.df$the.Mann.Whitney.Utest.combo <- as.character(new.df$the.Mann.Whitney.Utest.combo)
new.df$the.Mann.Whitney.Utest.combo <- as.numeric(new.df$the.Mann.Whitney.Utest.combo)
for.heat.map <- new.df[ which(  new.df$the.Mann.Whitney.Utest.combo <= 0.05),]
length(for.heat.map[,1]) # 86
genes.for.heat.map <- as.character(for.heat.map$the.gene)

low.CORT.M.subset <- low.CORT.M[,1]
names(low.CORT.M.subset)


# DECIDE IF WE WANT TO PLOT LOG2 VALUES OR NOT HERE ...
counts.for.heat.map <- cbind(low.CORT.F,low.CORT.M.subset,high.CORT)  
# counts.for.heat.map <-  2^(  counts.for.heat.map -0.25 )  # 
counts.for.heat.map <- counts.for.heat.map[genes.for.heat.map,] # subset by genes of interest

# rename
if(the.tissue == "cerebrum")
{colnames(counts.for.heat.map)[colnames(counts.for.heat.map) == 'low.CORT.M.subset'] <- 'MN21818Ce' }
if(the.tissue == "liver")
{colnames(counts.for.heat.map)[colnames(counts.for.heat.map) == 'low.CORT.M.subset'] <- 'MN21818Li' }
if(the.tissue == "kidney")
{colnames(counts.for.heat.map)[colnames(counts.for.heat.map) == 'low.CORT.M.subset'] <- 'MN21818Ki' }
if(the.tissue == "skeletalmuscle")
{colnames(counts.for.heat.map)[colnames(counts.for.heat.map) == 'low.CORT.M.subset'] <- 'MN21818Sm' }

length(colnames(counts.for.heat.map))



library("pheatmap")
mat  <- (counts.for.heat.map)[ genes.for.heat.map, ]
mat  <- mat - rowMeans(mat)

library(RColorBrewer)
cols <- blueColours <- rev(brewer.pal(9, "RdBu"))
# xx <- pheatmap(mat,annotation_col=anno, scale="row",color = cols)

# xx <- pheatmap(mat, scale="row",color = cols, clustering_distance_rows="euclidean") # with black border colour

xx <- pheatmap(mat, scale="row",color = cols, clustering_distance_rows="euclidean",border_color=NA)
# ?pheatmap






the.filename <-paste("CORT","_",the.tissue,"-heatmap.pdf",sep="")
# save_pheatmap_pdf <- function(x, filename, width=7, height=25) {
save_pheatmap_pdf <- function(x, filename, width=7, height=5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(xx, the.filename)








if(run.bootstrap == "yes")
{
  # assess bootstrap support
library(pvclust)
library(parallel)
# https://academic.oup.com/bioinformatics/article/22/12/1540/207339
# Compute p-value for hierarchical clustering
the.no.boostraps <- 100000
set.seed(123)
the.filename <-paste("CORT","_",the.tissue,"-bootstrap.pdf",sep="")
pdf(file = the.filename,   # The directory you want to save the file in
    width = 5, # The width of the plot in inches
    height = 4) # The height of the plot in inches
res.pv <- pvclust(mat, method.hclust = "complete",
        method.dist = "euclidean", nboot = the.no.boostraps , parallel=TRUE )
plot(res.pv, hang = -1, cex = 0.5)
pvrect(res.pv)
dev.off()
# Values on the dendrogram are AU p-values (Red, left), BP values (green, right), and clusterlabels (grey, bottom). Clusters with AU > = 95% are indicated by the rectangles and are considered to be strongly supported by dat
# pvclust provides two types of p-values: AU (Approximately Unbiased) p-value and BP (Bootstrap Probability) value. AU p-value, which is computed by multiscale bootstrap resampling, is a better approximation to unbiased p-value than BP value computed by normal bootstrap resampling.
}









# enable scientific notation
options(scipen=1)