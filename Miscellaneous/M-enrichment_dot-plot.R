# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Custom enrichment figure
#
# 2600520 - Inge Seim
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
rm(list=ls()) # reset R  

library(ggplot2)
library(openxlsx)



setwd("/Users/sciseim/Dropbox/seimlab\ at\ NNU/antechinus\ genome\ project/--MS--/REGs\ and\ PSGs/enrichment_dotplots/")




targetsheet <- 5 # GO BP





# Antechinus PSGs: interesting GO BP terms
the.enrichments <- read.xlsx("antechinuses PSGs Msigdb tables BGI-background with HPO.xlsx",sheet=targetsheet)
the.M <- 114 # no REGs or PSGs  
pathways.of.interest <- c("GO_REGULATION_OF_HUMORAL_IMMUNE_RESPONSE",
                          "GO_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_CIRCULATING_IMMUNOGLOBULIN",
                          "GO_LYMPHOCYTE_MEDIATED_IMMUNITY",
                          "GO_NEGATIVE_REGULATION_OF_IMMUNE_RESPONSE",
                          "GO_NEGATIVE_REGULATION_OF_INTERLEUKIN_17_PRODUCTION",
                          "GO_NEGATIVE_REGULATION_OF_ACUTE_INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS",
                          "GO_RESPONSE_TO_INTERLEUKIN_3",
                          
                          "GO_GLYCOLIPID_TRANSPORT",
                          "GO_NEUTRAL_LIPID_CATABOLIC_PROCESS",
                          
                          "GO_BEHAVIORAL_RESPONSE_TO_PAIN",
                          "GO_RESPONSE_TO_PAIN",
                          
                          
                          "GO_COMPLEMENT_ACTIVATION",
                          "GO_NEGATIVE_REGULATION_OF_COMPLEMENT_ACTIVATION",
                          "GO_COMPLEMENT_ACTIVATION_LECTIN_PATHWAY"
)
























# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
the.enrichments.SUBSET <- the.enrichments[ which(the.enrichments$Pathway.Name %in% pathways.of.interest) ,]
names(the.enrichments.SUBSET)
# calculate enrichment. see http://metascape.org/blog/?p=122
# There are N total number of genes in our study pool (this is also known as the “background” gene list, defaults to all genes in the genome). 
# A given pathway of interest consists of k gene members. 
#Our input gene list consists of M genes, among which n are found to fall into the same given pathway.
the.N <- 6255  # 
# the.k <- # pathways in our background
# the.M <- 455 # no REGs or PSGs  
# the.n <- # input genes in pathway  
# enrichment <- nN/kM
the.enrichments.SUBSET$enrichments <- (the.enrichments.SUBSET$Observed.Number * the.N) / (the.enrichments.SUBSET$Pathway.Size * the.M)  

enrichments <- the.enrichments.SUBSET$enrichments
pvalues <- the.enrichments.SUBSET$Bootstrap.P.Value
pathways <- the.enrichments.SUBSET$Pathway.Name
nogenes <- the.enrichments.SUBSET$Observed.Number




the.data <- as.data.frame( cbind(enrichments,pvalues,pathways,nogenes) )
the.data$enrichments  <- as.numeric(the.data$enrichments)
the.data$log.enrichment <- log ( the.data$enrichments )
class(the.data)



# pathways <- gsub("GO_","",pathways) # get rid of prefixes
# should create new y-axis labels
# pathways <- paste(gsub("GO_","",pathways)," (",nogenes,"/",the.M,")",sep="")
pathways <- paste(gsub("GO_","",pathways)," (",nogenes,")",sep="")
the.data$pathways <- pathways
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# sort by general function









# Antechinus PSGs
target.order <- the.data$pathways # keep same
target.order <- c(
  "NEUTRAL_LIPID_CATABOLIC_PROCESS (2)",                                         
  "GLYCOLIPID_TRANSPORT (1)",                                                    
  
  "BEHAVIORAL_RESPONSE_TO_PAIN (2)",                                             
  "RESPONSE_TO_PAIN (2)",                                                        
  
  "COMPLEMENT_ACTIVATION (3)",                                                   
  "NEGATIVE_REGULATION_OF_COMPLEMENT_ACTIVATION (1)",                            
  "COMPLEMENT_ACTIVATION_LECTIN_PATHWAY (2)",                                    
  
  "NEGATIVE_REGULATION_OF_IMMUNE_RESPONSE (3)",                                  
  "HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_CIRCULATING_IMMUNOGLOBULIN (3)",          
  "REGULATION_OF_HUMORAL_IMMUNE_RESPONSE (2)", 
  "LYMPHOCYTE_MEDIATED_IMMUNITY (4)",                                            
  "NEGATIVE_REGULATION_OF_INTERLEUKIN_17_PRODUCTION (1)",                        
  "NEGATIVE_REGULATION_OF_ACUTE_INFLAMMATORY_RESPONSE_TO_ANTIGENIC_STIMULUS (1)",
  "RESPONSE_TO_INTERLEUKIN_3 (1)"
)























# can skip this and force order in the actual ggplot by assigning a level!
 row.names(the.data) <- the.data$pathways
 the.data <- the.data[target.order,]
# and to make R respect the order
 the.data$pathways <- factor(the.data$pathways , levels = the.data$pathways [order(the.data$pathways )])

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
  
 
 enrichment.data.for.plot <- the.data$log.enrichment
 enrichment.data.for.plot <- the.data$enrichments
 enrichment.data.for.plot.no.log <- the.data$enrichments
 
 the.data$pvalues  <- as.numeric( as.character(the.data$pvalues) )
 class(the.data$pvalues) # numeric
 
 
 
 
# testing 
library(viridis)
ggplot(data=the.data, aes(x=pvalues, y=pathways, color=enrichment.data.for.plot)) + geom_point(alpha = 0.8,size=5) + theme_classic()  + labs(x = "P-value", y = "",colour = "Enrichment\n") + scale_color_viridis(option = "D") + theme(text = element_text(size=16)) + scale_x_continuous(name = "P-value", limits = c(0, 0.05))
# View(the.data)

# try to increase the font sizes of the axes
# theme(text = element_text(size=20)

library("ggplot2")
S1 <- ggplot(data=the.data, aes(x=pvalues, y=factor(pathways,level=target.order), color=enrichment.data.for.plot)) + geom_point(alpha = 0.8,size=5) + theme_classic()  + labs(x = "P-value", y = "",colour = "Enrichment\n") + theme(text = element_text(size=16)) + scale_x_continuous(name = "P-value", limits = c(0, 0.05))
# + theme(axis.text.x=element_text(size=rel(0.5)))






library(viridis)
# ?viridis
# see https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html#gallery
S1 + scale_color_viridis(option = "D")
# S1 + scale_color_viridis(option = "magma")

# if we want to resize by no genes
#@ S1 <- ggplot(data=the.data, aes(x=pvalues, y=pathways, color=enrichment.data.for.plot)) + geom_point(alpha = 0.8,size=nogenes) + theme_classic()  + labs(x = "P-value", y = "",colour = "Enrichment\n")
#@ S1 + scale_color_viridis(option = "D")

the.enrichment.plot <- S1 + scale_color_viridis(option = "D")
pdf("-enrichment_plot.PDF",width=13)
# 20 for PSG, 10 for REGs

the.enrichment.plot
dev.off()

