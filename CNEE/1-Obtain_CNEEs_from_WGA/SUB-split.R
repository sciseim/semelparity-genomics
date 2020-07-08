big.FASTA.split <- strsplit(big.FASTA,">human")
# head(big.FASTA.split)
big.FASTA.split.character <- unlist(big.FASTA.split)
# class(big.FASTA.split.character)

# freem up memory
remove(big.FASTA)
remove(big.FASTA.split)


# bad R code since so much in memory, but works 
big.FASTA.split.character <- big.FASTA.split.character[2:length(big.FASTA.split.character)] # first element is blank. Remove it!



# system("mkdir -p ./split")

# write each to new file
# loop through string

for (CNEE in big.FASTA.split.character) 
{
  CNEE <- paste(">human", CNEE,sep="")   # add human fasta header back
  CNEE.counter = CNEE.counter + 1
  the.filename  <- paste("../split/CNEE",CNEE.counter,".fa",sep="")
  cat(CNEE,file=the.filename)
}

remove(big.FASTA.split.character)