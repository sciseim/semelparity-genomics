# install.packages(c("stringr"))
# install.packages(c("seqinr"))
# install.packages(c("dplyr"))
library("stringr")
library("seqinr")
library("dplyr")


# run for each chromosome 
setwd("~/Downloads/PhyloAcc_prep/split")
system("mkdir -p ../no_missing_species")





# START BIG LOOP
files <- list.files(pattern = "\\.fa$")
for(i in files){
  # fasta.name <- "CNEE17.fa" # has 19 
  fasta.name <- i
  fastafile <- read.fasta(file = paste("./",fasta.name, sep=""),seqtype = "DNA",as.string = TRUE, set.attributes = FALSE) 
  
  names(fastafile)
  
  # get rid of suffixes ... i.e. only keep human
  names(fastafile) <- sub("\\..*", "", names(fastafile))
  
  
  
  
  desired.number.of.species <- 19
  if(length(names(fastafile)) == desired.number.of.species )
  {
  cat("PASS!")  
  the.CNEE.name <- gsub("*.fa","",fasta.name)  


# write out if CNEE has all 19 species
# remove any =
  CNEE.outputsequence <- fastafile
  names(CNEE.outputsequence)
  length(fastafile) # 2
  CNEE.outputsequence[19] <- gsub("[=]","",CNEE.outputsequence[19])    
  length(CNEE.outputsequence) # 19


    
# only output if length cutoff passes
# system("mkdir -p ../no_missing_species")
CNEE.length.cutoff <- 30
if(    nchar(CNEE.outputsequence[1])  >= CNEE.length.cutoff )
{  
file.out <- paste("../no_missing_species/",fasta.name,sep="")
write.fasta(CNEE.outputsequence, names(CNEE.outputsequence), file.out, open = "w", nbchar = 60, as.string = FALSE)
# not appending here, so 'w'
}  


# the below will output files used to make 100% sure the output and bed files are correct!    
# write to FASTA *per* species
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (species.number in 1:desired.number.of.species) {
  outputspecies <- paste("species_",species.number,sep="")
  outputsequence <- fastafile[species.number]
  
  # remove any =
  outputsequence <- gsub("[=]","",outputsequence)  
  
  
  names(outputsequence) <- the.CNEE.name
  file.out <- paste("../",outputspecies,".fa",sep="")
  write.fasta(outputsequence, names(outputsequence), file.out, open = "a", nbchar = 60, as.string = FALSE)
  # append here since we will add all CNEEs as separate entries for each species
}  
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
}
   
  
  
}
# END OF BIG LOOP
   
  





