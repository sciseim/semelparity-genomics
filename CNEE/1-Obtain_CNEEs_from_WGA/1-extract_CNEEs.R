rm(list = ls())
library(readr)

# start count here so that we can parse through multiple MAFs...
# cannot load more than 2Gb into RAM in R currently, so this is the only way to ensure that we obtain CNEEs in order when different chromosomes are parsed
CNEE.counter = 0

setwd("~/Downloads/PhyloAcc_prep/input/")
system("mkdir -p ../split")


# cat NC_000001.11.maf | awk '/^s/{print ">" $2 "\n" $7}' >NC_000001.11.maf.bigFA



# run for each chromosome
# big.FASTA <- read_file("NC_000024.10.bigFA") # chrY ... 2,900 CNEEs ... 4.5Mb
# big.FASTA <- read_file("NC_000001.11.maf.bigFA") # chr1 ... 109,877 CNEEs ... 204.8 Mb

# before filtering, we have x CNEEs 



files <- list.files(pattern = "\\.fa$")
for(i in files)
{
  # cat("add to empty DF\n")
  fasta.name <- i
cat(fasta.name,"\n")
big.FASTA <- read_file(fasta.name) # load FASTA
#
 source("../SUB-split.R")

}

