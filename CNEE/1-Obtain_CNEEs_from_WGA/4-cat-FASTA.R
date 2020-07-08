rm(list=ls()) 
library("seqinr")

# setwd("~/Downloads/obtain CNEEs/split/no_missing_species/")
setwd("~/Downloads/PhyloAcc_prep/no_missing_species/")


system("rm -rf ./tmp")
system("mkdir ./tmp")
system("rm ./PhyloAccInout/input.fa") # delete the output file



# load FASTA files

files <- list.files(pattern = "\\.fa$")

# load random FASTA and grab headers (all the same)
fasta.name <- files[1]
fastafile <- read.fasta(file = paste("",fasta.name, sep=""),seqtype = "DNA",as.string = TRUE, set.attributes = FALSE) 
TARGETSPECIES <- names(fastafile)
# for each species... pull out its sequence from each FASTA




class(TARGETSPECIES) # character
for (the.targetspecies in TARGETSPECIES)
{ 
  # START OF TARGET SPECIES LOOP
  # chromosomes[i]
  print(the.targetspecies)
  file.name <- paste(the.targetspecies,".fa",sep="")
  
  # pull out from each FASTA
  for(i in files)
  {
    # cat("add to empty DF\n")
    
    # i <- "VCE10084.fa"  
    # i <- "VCE10138.fa"  
    fasta.name <- i  
    fastafile <- read.fasta(file = paste("",fasta.name, sep=""),seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)   
    targetseq.fasta <- fastafile[which(names(fastafile) == the.targetspecies)  ]
    targetseq.fasta <- targetseq.fasta[[1]]  # convert to string
    
    # save output. String only... load back in and add header after all sequences added  
    cat(targetseq.fasta,file=paste("./tmp/",the.targetspecies,".txt",sep=""),sep="",append=TRUE)
    
  } # END OF FASTA LOOP
} # END OF SPECIES LOOP



# then prepare final file
for (the.targetspecies in TARGETSPECIES)
{ 
  # START OF TARGET SPECIES LOOP
  # chromosomes[i]
  print(the.targetspecies)
  
  fileName <- paste("./tmp/",the.targetspecies,".txt",sep="")
  cat.ed.FASTA.string <- paste(readLines(fileName), collapse=" ")
  # class(cat.ed.FASTA.string)
  
  # 150819 - replace N with -
  cat.ed.FASTA.string <- gsub("N","-",cat.ed.FASTA.string)
  cat.ed.FASTA.string <- gsub("n","-",cat.ed.FASTA.string)
  
  
  fasta.name <- the.targetspecies
  file.out <- paste("./PhyloAccInout/input.fa",sep="")
  
  write.fasta(cat.ed.FASTA.string, fasta.name, file.out, open = "a", nbchar = 60, as.string = FALSE)
  # "a" here -- want to append
  
}






