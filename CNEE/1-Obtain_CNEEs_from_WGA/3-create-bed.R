rm(list=ls()) 


# setwd("/Volumes/KAMDISK/PhyloAcc/My_PhyloAcc/Step7/X-post-mafft-filtering/CNEEs_mafft_trimAl/subset/CNEE_no_missing_species")
# setwd("~/Downloads/obtain CNEEs/split/no_missing_species/")
setwd("~/Downloads/PhyloAcc_prep/no_missing_species/")


system("rm ./PhyloAccInout/input.bed") # delete the output file


# PhyloAcc input
# (2) a multiple alignment files concatenating sequences of all input conserved elements in FASTA format, and 
#
# (3) a bed file with the position of each individual element in the coordinate of concatenated alignment file (0-based).

# for each FASTA file ...
# add to single CNEE per SPECIES


# load FASTA files


files <- list.files(pattern = "\\.fa$")



# read in each
remove(CNEEINFO)
CNEEINFO = NULL # create empty df
for(i in files)
{
  # cat("add to empty DF\n")
  fasta.name <- i
  
  # testing
  # fasta.name <- "VCE10084.fa"
  # fasta.name <- "VCE10138.fa"
  
  fastafile <- read.fasta(file = paste("",fasta.name, sep=""),seqtype = "DNA",as.string = TRUE, set.attributes = FALSE) 
  CNEE.name <- gsub(".fa","",fasta.name)
  hg38.CNEE.seq <- (fastafile[[which(names(fastafile) == "human")]])
  CNEE.length <- nchar(hg38.CNEE.seq)
  CNEE.filename <- fasta.name
  
  
  temp.df <- c(CNEE.name,CNEE.filename,CNEE.length)
  CNEEINFO <- rbind(CNEEINFO,temp.df)  
}



# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# next obtain coordinates and write a FASTA file + bed

# SAVE OUTPUT 
# getwd() # "/Step7/X-post-mafft-filtering/CNEEs_mafft_trimAl"
system("rm -rf ./PhyloAccInout")
system("mkdir -p ./PhyloAccInout")

temp <- as.data.frame(CNEEINFO)
names(temp) <- c("CNEE","filename","CNEE.length")
head(temp)

# save first CNEE first

# SEG_FILE: the path of bed file for genomic regions. At least 3 columns in the bed file. You could have more columns in this file but the program will only read in the first 3 columns. The first column is element name or ID (different from usual bed file), the second and third columns are the starting and ending positions of each element (in the coordinate of the whole multiple alignment). The program assumes that the alignment file concatenates all the elements together and will only use the second and third columns in the bed file. If concatenating multiple chromosomes, the coordinate of elements on the current chromosome should not start from zero but should add to the previous chromosome. The program will internally generated a No. for each element which is the order in the input bed file, and it will use No. in the outputs and plot functions.
CNEE.name.column <- as.character(temp[(1),1])
CNEE.name.column <- gsub("VCE","",CNEE.name.column) # convert VCE10084 to 10084

first.bed.data <- paste(CNEE.name.column,1,as.numeric(as.character(temp[(1),3])),sep="\t")
bed.file.name <- paste("./PhyloAccInout","/input.bed",sep="")
write.table(first.bed.data, file =bed.file.name, sep = "\t",
            row.names = FALSE,quote = FALSE,append = FALSE,col.names = FALSE)

# loop through all after the first sequence
# new.temp.df.row= NULL
for(j in 2:nrow(temp)) {
  # for each row from row 2 to the end
  # j <- 2
  the.row <- temp[j,]
  
  # start.this.CNEE <- as.numeric(as.character((temp[(j),3])))
  size.of.this.CNEE <- as.numeric(as.character((the.row[,3]))) # save as above
  
  # need to sum up all CNEEs below
  head(temp)
  
  end.of.CNEE.below <- sum(as.numeric(as.character((temp[1:(j-1),3]))))
  
  # end.of.CNEE.below <- as.numeric(as.character((temp[(j-1),3])))
  the.row$bedFROM <- end.of.CNEE.below+1
  the.row$bedTO <- the.row$bedFROM+size.of.this.CNEE-1
  
  CNEE.name.column <- as.character(the.row$CNEE)
  CNEE.name.column <- gsub("VCE","",CNEE.name.column) # convert VCE10084 to 10084
  
  the.bed.data <- paste(CNEE.name.column,as.numeric(as.character(the.row[(1),4])),as.numeric(as.character(the.row[(1),5])),sep="\t")
  
  
  bed.file.name <- paste("./PhyloAccInout","/input.bed",sep="")
  write.table(the.bed.data, file =bed.file.name, sep = "\t",
              row.names = FALSE,quote = FALSE,append = TRUE,col.names = FALSE)
  # and append!
}

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
