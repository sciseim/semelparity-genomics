rm(list=ls()) 

setwd("~/Downloads/pseudopipe results/")

fileprefix <- "G.agilis.fushu" # A.flav.daiqu or G.agilis.fushu
thefile <-paste(fileprefix,".gene.transcript.txt",sep="")
pseudopipe.output <- read.delim2(thefile,header=FALSE)
# and the output
outputfile <- paste(fileprefix,"-pseudopipe.output.class.shifts_et_stops.idcutoff.txt",sep="")

head(pseudopipe.output)

#chr start end strand query frac ins del shift stop expect ident polya type
names(pseudopipe.output) <- c("chr", "start","end", "strand" ,"query" ,"frac", "ins", "del", "shift" ,"stop" ,"expect" ,"ident" ,"polya", "type")
# https://faq.gersteinlab.org/2019/06/03/terms-in-pseudopipe-output-etc/
# frac = fraction of parent gene that matches the pseudogene
# ins = number of insertions in the pseudogene compared to parent sequence
# del = number of deletions in the pseudogene compared to parent sequence
# shift = number of frame shifts in the pseudogene compared to parent sequence
# stop = number of stop codons in the pseudogene compared to parent sequence
# polya = flag indicating the presence or absence of a polyA tail **many true PSSDs will contain this!***

# PseudoPipe aligns the query proteins against ***regions of the genome that are not already annotated as genic regions***, and then classifies the resulting protein sequences as a variety of possible pseudogenes, including “retrotransposed pseudogenes” (“PSSD,” which lack intronic sequence), “duplicated pseudogenes” (“DUP,” which may contain intronic sequence) and “pseudogenic fragments” (“FRAG,” which are too fragmented to determine their source).

# Pseudogenic fragments (FRAGs) = are fragments that have high-sequence similarity to known proteins, but are too decayed to be reliably assessed as processed or duplicated. 


# PseudoPipe35 and Retrofinder36 pipelines (details in Methods). These computational pseudogene predictions provide ***hints*** to manual annotators during the first-pass of annotation and identify potential missing features, flagging them for manual re-investigation (Fig. 1).
# from https://www.nature.com/articles/nature28172.pdf?proof=trueMay.


# Unitary pseudogenes are class of unprocessed pseudogenes without functioning (i.e. genic) counterparts in the genome.
# Numerically, they constitute only a small fraction of tens of thousands of annotated pseudogenes in the human genome.
# 



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# columns currently factors -- change to correct format
head(pseudopipe.output)
pseudopipe.output$start <- as.numeric(as.character(pseudopipe.output$start))
pseudopipe.output$end <- as.numeric(as.character(pseudopipe.output$end))
pseudopipe.output$frac <- as.numeric(as.character(pseudopipe.output$frac))
pseudopipe.output$ins <- as.numeric(as.character(pseudopipe.output$ins))
pseudopipe.output$del <- as.numeric(as.character(pseudopipe.output$del))
pseudopipe.output$shift <- as.numeric(as.character(pseudopipe.output$shift))
pseudopipe.output$stop <- as.numeric(as.character(pseudopipe.output$stop))
pseudopipe.output$expect <- as.numeric(as.character(pseudopipe.output$expect))
pseudopipe.output$ident <- as.numeric(as.character(pseudopipe.output$ident))
pseudopipe.output$polya <- as.numeric(as.character(pseudopipe.output$polya))
pseudopipe.output$chr <- as.character(pseudopipe.output$chr)
pseudopipe.output$strand <- as.character(pseudopipe.output$strand)
pseudopipe.output$query <- as.character(pseudopipe.output$query)
pseudopipe.output$type <- as.character(pseudopipe.output$type)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# decide on cut-off
# PGD: a pangolin genome hub for the research community
# https://academic.oup.com/database/article/doi/10.1093/database/baw063/2630389
# (1E−10 e-value, 70% parent gene coverage, 40% gene identity)
# -10 possibly too high here since Tassie devil is the input species (above is vs. self in pangolin)

table(pseudopipe.output$type)
# DUP FRAG PSSD 
# 533 3783 1898
length (  which(pseudopipe.output$type == "DUP")  ) # duplicated
length (  which(pseudopipe.output$type == "PSSD")  ) # processed pseudogene
length (  which(pseudopipe.output$type == "FRAG")  ) # fragments

# keep hits that contain at least 70% of parent transcript
pseudopipe.output <- pseudopipe.output[  (which(pseudopipe.output$frac >= 0.50)) , ] # 2,540
length(pseudopipe.output$chr) # 2,540

# only keep fragments
pseudopipe.output.FRAG <- pseudopipe.output[  (which(pseudopipe.output$type == "FRAG")) , ]
length(pseudopipe.output.FRAG$chr) # 302

# keep fragments *and* dup
pseudopipe.output.class <- pseudopipe.output[  (which(pseudopipe.output$type == "FRAG" | pseudopipe.output$type == "DUP")) , ]
length(pseudopipe.output.class$chr) # 642





# if a polyA, likely a processed pseudogene, so probably a genuine unitary pseudogene (issue for FRAG class and, of course, single-exon DUPL class)!
pseudopipe.output.class <- pseudopipe.output.class[  (which(pseudopipe.output.class$polya == 0)) , ]
length(pseudopipe.output.class$chr) # 454


length(   (which(pseudopipe.output.class$stop >= 1)) ) # 394
length(   (which(pseudopipe.output.class$shift >= 1)) ) # 392

# let us consider candidate pseudogenised genes with >=1 stop codon and/or >=1 frameshift
pseudopipe.output.class.shifts_et_stops <- pseudopipe.output.class[  (which(pseudopipe.output.class$stop >= 1 | pseudopipe.output.class$shift >= 1)) , ] 
length(pseudopipe.output.class.shifts_et_stops$chr) # 416


# identity above
pseudopipe.output.class.shifts_et_stops.idcutoff <- pseudopipe.output.class.shifts_et_stops[  (which(pseudopipe.output.class.shifts_et_stops$ident >= 0.30)) , ] # 
length(pseudopipe.output.class.shifts_et_stops.idcutoff$chr) # 342 at 40%
length(pseudopipe.output.class.shifts_et_stops.idcutoff$chr) # 410 at 30%

# 


pseudopipe.output.class.shifts_et_stops.idcutoff$ident 
pseudopipe.output.class.shifts_et_stops.idcutoff




























# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# TasDev gff
# load devil information from gff
the.devil.wears.gffs <- read.delim2("PARSED_devil_gff-in_loop.tsv",header=FALSE)
length(the.devil.wears.gffs[,1]) # number of genes (lines)

# e.g. SEPTIN7	XM_023498535.2	XM_023498535.2	XP_023354303.2	septin-7 isoform X3
names(the.devil.wears.gffs) <- c("gene","RefSeq.RNA.loop","RefSeq.RNA","RefSeq.protein","description")
# can drop the RNA.loop (just a control for 'X-parse-devil-gff3-v2.R')
the.devil.wears.gffs <- the.devil.wears.gffs[c("gene","RefSeq.RNA","RefSeq.protein","description")]


length(unique(the.devil.wears.gffs$gene)) # unique genes (have 23,124 protein-coding genes and a total of 45k + transcripts)
length(unique(the.devil.wears.gffs$RefSeq.RNA)) # 



head(the.devil.wears.gffs)

pseudopipe.output.class.shifts_et_stops.idcutoff

head(pseudopipe.output.class.shifts_et_stops.idcutoff)
# look-up
merge(merge(pseudopipe.output.class.shifts_et_stops.idcutoff, the.devil.wears.gffs, by.x = "query",by.y="RefSeq.protein")[1:3])


pseudopipe.temp <- pseudopipe.output.class.shifts_et_stops.idcutoff
pseudopipe.temp$RefSeq.protein <- pseudopipe.output.class.shifts_et_stops.idcutoff$query
# Use plyr to rename column and to mergge (plyr join has more function)
# http://www.inside-r.org/packages/cran/plyr/docs/join
library(plyr) #  from CRAN
# can now join!
pseudopipe.final <- join(pseudopipe.temp, the.devil.wears.gffs, type = "left", match = "first")
head(pseudopipe.final)
length(pseudopipe.final[,1])
# View(pseudopipe.final)
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# save output ... move to bottom once devil GFF3 has been parsed
# head(pseudopipe.final)
# write.table(pseudopipe.final, outputfile, sep="\t",row.names = FALSE,append=FALSE,col.names = TRUE,quote = FALSE)
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@















# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# only keep pseudogenes that do NOT have a homolog in the gene annotation (BGI EVM)
#
#
# even though we *did* input a gff to ignore annotated genes, the pipeline will may identify pseudogenes of annotated genes in our query species => filter again to keep genes NOT annotated for our species of interest!


# maybe load gene models so that we only validate candidate pseudogenes *without* a gene model in our semelparous species
# the.evm.gene.models <- read.delim2("A.flav-the.annotations.txt",header=FALSE)
the.evm.gene.models <- read.delim2(paste(fileprefix,"-the.annotations.txt",sep=""),header=FALSE)
length(the.evm.gene.models[,1]) # number of genes (lines) - 24,707
#
#
not.in.EVM.genes <- setdiff(pseudopipe.final$gene, the.evm.gene.models[,2])
pseudopipe.final <- pseudopipe.final[   which(pseudopipe.final$gene %in% not.in.EVM.genes) , ] 
# View(pseudopipe.final)
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@










# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# save output ... move to bottom once devil GFF3 has been parsed
#
head(pseudopipe.final)
write.table(pseudopipe.final, outputfile, sep="\t",row.names = FALSE,append=FALSE,col.names = TRUE,quote = FALSE)
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@























