


#  N for unknown nucleic acid residue or X for unknown amino acid residue
# replace with a gap AKA '-'
replace.NorX.with.gap <- "yes"
# yes 
# gapped with -/dashes for alignments
# e.g. ATG---CATCATCAT
# e.g. M-HHH
#
# no 
# gapped with n or X for BLAST db creation
# same as default gffread output (expect lowercase) -- same output as Tasmanian tiger pseudo-sequences
# e.g. ATGnnnCATCATCAT
# e.g. MxHHH


the.samples <- c("BD-17-5A","42595A","46020A","AA100A")
# the.sample <- "BD-17-5A"
# Murexia sp. = 42595A
# Murexia melanurus = 46020A
# Antechinus arktos = AA100A
# Antechinus argentus = BD-17-5A

# RUN FOR EACH SAMPLE
for (the.sample in the.samples) {
  




directory.path <- paste("~/Downloads/testar/",the.sample,sep="")
setwd(directory.path)
# system("mkdir -p READY_FOR_ANALYSIS")




if(the.sample == "42595A")
{
  transcript.output <- "Murexia_sp-TasDevPseudo-RNA.fa"
  CDS.output <- "Murexia_sp-TasDevPseudo-CDS.fa"
  AA.output <- "Murexia_sp-TasDevPseudo-PROTEIN.fa"
  species.header <- "M_sp"
}
if(the.sample == "46020A")
{
  transcript.output <- "Murexia_melanurus.-TasDevPseudo-RNA.fa"
  CDS.output <- "Murexia_melanurus-TasDevPseudo-CDS.fa"
  AA.output <- "Murexia_melanurus-TasDevPseudo-PROTEIN.fa"
  species.header <- "M_melanurus"
}
if(the.sample == "BD-17-5A")
{
  transcript.output <- "Antechinus_argentus-TasDevPseudo-RNA.fa"
  CDS.output <- "Antechinus_argentus-TasDevPseudo-CDS.fa"
  AA.output <- "Antechinus_argentus-TasDevPseudo-PROTEIN.fa"
  species.header <- "A_argentus"
}
if(the.sample == "AA100A")
{
  transcript.output <- "Antechinus_arktos-TasDevPseudo-RNA.fa"
  CDS.output <- "Antechinus_arktos-TasDevPseudo-CDS.fa"
  AA.output <- "Antechinus_arktos-TasDevPseudo-PROTEIN.fa"
  species.header <- "A_arktos"
}








# on server
# system("cat ./CDS/*.fa >CDS.mfa")
# system("cat ./AA/*.fa >AA.mfa")
# system("cat ./transcript/*.fa >transcripts.mfa")
# cat ./CDS/*.fa >CDS.mfa && cat ./AA/*.fa >AA.mfa && cat ./transcript/*.fa >transcripts.mfa


# install.packages("seqinr")
library("seqinr")
the.transcripts <- read.fasta(file = "transcripts.mfa")
the.CDS <- read.fasta(file = "CDS.mfa")
the.AA <- read.fasta(file = "AA.mfa")

#how many fasta sequences
length(the.transcripts) # 24,041
length(the.CDS) # 22,372
length(the.AA) # 22,372




# #####################
# CDS
# START OF BIG LOOP
length(the.CDS) 
for (i in 1:(length(the.CDS))) 
{
  # i <- 1 # dummy
  targetseq <- the.CDS[[i]] # # Put the sequence in a vector
  GC(targetseq)
  length(targetseq) # 1699 bp
  length(which(targetseq %in% "n"))
  
  if(length(which(targetseq %in% "n")) != length(targetseq))
  {"output!"
    class(targetseq) 
    the.name <- getName(targetseq)
    the.name <- paste(the.name,"_",species.header,sep="") 

    
    if(replace.NorX.with.gap == "yes")
    {
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    # 15.05.19
    # replace 'n' with dash/-
    class(targetseq) # SeqFastadna
#    which(targetseq == "n") 
#    length(which(targetseq == "n")) # could be >1!
    #
    targetseq <- gsub('n', '-', targetseq)
    targetseq <- gsub('>', '-', targetseq)
    #    length(which(targetseq == "n")) # could be ZERO
#    length(which(targetseq == "-")) # could be same as before 'gsub'
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    }
    
    
    # write to FASTA (append)
    write.fasta(targetseq, the.name, paste(CDS.output,sep=""), open = "a", nbchar = 60, as.string = FALSE)
    # a=append https://www.rdocumentation.org/packages/seqinr/versions/3.4-5/topics/write.fasta
  }
} # END OF BIG LOOP
# #####################

# #####################
# transcripts
# START OF BIG LOOP
length(the.transcripts) 
for (i in 1:(length(the.transcripts))) 
{
  # i <- 1 # dummy
  targetseq <- the.transcripts[[i]] # # Put the sequence in a vector
  GC(targetseq)
  length(targetseq) # 1699 bp
  length(which(targetseq %in% "n"))
  
  if(length(which(targetseq %in% "n")) != length(targetseq))
  {"output!"
    class(targetseq) 
    the.name <- getName(targetseq)
    the.name <- paste(the.name,"_",species.header,sep="") 
    
    
    if(replace.NorX.with.gap == "yes")
    {
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    # 15.05.19
    # replace 'n' with dash/-
    class(targetseq) # SeqFastadna
#    which(targetseq == "n") 
#    length(which(targetseq == "n")) # could be >1!
    #
    targetseq <- gsub('n', '-', targetseq)
    targetseq <- gsub('>', '-', targetseq)
    length(which(targetseq == ">")) # could be >1!
#    length(which(targetseq == "n")) # could be ZERO
#    length(which(targetseq == "-")) # could be same as before 'gsub'
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    }

        
    # write to FASTA (append)
    write.fasta(targetseq, the.name, paste(transcript.output,sep=""), open = "a", nbchar = 60, as.string = FALSE)
    # a=append https://www.rdocumentation.org/packages/seqinr/versions/3.4-5/topics/write.fasta
  }
} # END OF BIG LOOP
# #####################

# #####################
# AA
# START OF BIG LOOP
length(the.AA) 
for (i in 1:(length(the.AA))) 
{
  # i <- 1 # dummy
  targetseq <- the.AA[[i]] # # Put the sequence in a vector
  GC(targetseq)
  length(targetseq) # 1699 bp
  length(which(targetseq %in% "n"))
  
  if(length(which(targetseq %in% "n")) != length(targetseq))
  {"output!"
    class(targetseq) 
    the.name <- getName(targetseq)
    the.name <- paste(the.name,"_",species.header,sep="") 
  
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    
    if(replace.NorX.with.gap == "yes")
    {
    # 15.05.19
    # replace 'n' with dash/-
    class(targetseq) # SeqFastadna
#    which(targetseq == "X") 
#    length(which(targetseq == "X")) # could be >1!
    #
    targetseq <- gsub('x', '-', targetseq)
    #    length(which(targetseq == "X")) # could be ZERO
#    length(which(targetseq == "-")) # could be same as before 'gsub'
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
    }
    
    
    # write to FASTA (append)
    write.fasta(targetseq, the.name, paste(AA.output,sep=""), open = "a", nbchar = 60, as.string = FALSE)
    # a=append https://www.rdocumentation.org/packages/seqinr/versions/3.4-5/topics/write.fasta
  }
} # END OF BIG LOOP
# #####################



} # END OF FOR EACH SAMPLE LOOP
