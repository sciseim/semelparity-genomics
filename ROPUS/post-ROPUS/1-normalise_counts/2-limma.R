
# setwd("../")
curated.dataset <- readRDS("curated_dataset.rds")
counts.LOG <- readRDS("normalized_counts.rds")
counts <- 2^(   counts.LOG )-0.25
# counts.LOG <- log2(counts + 0.25)


# create output dir
system("mkdir -p ./output")



# the.tissues <- unique(curated.dataset$tissue)
the.tissues <- c("cerebrum","gastrointestine","kidney","liver","skeletalmuscle","spleen","stomach","sternalgland","adrenalgland")


target.tissue <- i

for (i in the.tissues)
{
target.tissue <- i
#  target.tissue <- "liver"
try(source("SUB.limma.R"))
}


