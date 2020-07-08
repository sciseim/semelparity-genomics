#!/usr/bin/env Rscript 

# to run: from Rscript ./one_GTF_per_gene.Rscript $REFERENCEGTFPATH
args <- commandArgs()
print(args)


# should be separate files (GTF and later CDS/AA/transcript)
# evm.model.ctg2_pilon_pilon.3
# evm.model.ctg2_pilon_pilon.30
# evm.model.ctg2_pilon_pilon.300


# you will need data.table
# we will use fread to quickly load the very large GTF into RAM
# install.packages(data.table)

# setwd("~/Downloads/")
the.gtf.path <- name <- args[6]

# the.gtf <- read.delim2(the.gtf.path,header=FALSE,quote = "#")
# too slow

# output directory
system("mkdir -p singleGTFs")


library(data.table)
the.gtf <- fread(the.gtf.path,header=FALSE)
# the.gtf <- fread("yellowfoot_genome.FINAL.from_draft.masked.gtf",header=FALSE)
# the.gtf <- fread("dummy.gtf",header=FALSE)
# it handles doublequote automatically e.g. #!

# colnames(the.gtf)
# the.gtf[1:1,7:9]
# head(the.gtf)
# class(the.gtf)


# dummy until I am 100% sure it is working
# the.gtf <- the.gtf[1:10000,] # Testing: 347 genes
# tail(the.gtf)
# ENSSHAG00000016243 ... should be at least three genes now... is 5

# chr1    EVM    CDS    103067    103381    .    -    0    transcript_id "evm.model.ctg48_pilon_pilon.1";
# Antechinus
# chr: chr1
# source: EVM
# feature: CDS
# start: 103067
# end: 103381
# score:  .
# strand: -
# frame: 0
# attribute: transcript_id
# transcript_id "evm.model.ctg48_pilon_pilon.1";

colnames(the.gtf) <- c("chr","source","feature","start","end","score","strand","frame","attribute")
# the.gtf[1:1,1]
# class(the.gtf$attribute)
the.gtf$attribute <- as.character(the.gtf$attribute)


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# BIG LOOP PER CHR/SCAFFOLD/CONTIG HERE
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# NEED TO SPLIT INTO DFs PER SCAFFOLD
# the.gtf[1:1,1]
# chr
# 1: GL849905.1
the.chr <- the.gtf$chr
the.chr <- as.character(the.chr)
# length(the.chr) # 1000 (dummy)
unique.chr <-  the.chr[!duplicated(the.chr)]
length(unique.chr) # 43 scaffolds

# i <- 8 # dummy A. flavipes
for (i in 1:(length(unique.chr)))
# START OF BIG LOOP
{
    scaffold.of.interest <- unique.chr[i]
    # GL834541.1
    keep.these <- which(the.gtf$chr %in% scaffold.of.interest)
    the.gtf.subset.TMP <- the.gtf[keep.these,]
    
    
    
    
    
    
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # here comes the per scaffold work
    # the.gtf.subset.TMP$attribute
    
    # devil
    # [1000] gene_id "ENSSHAG00000018951"; gene_version "1"; transcript_id "ENSSHAT00000022567"; transcript_version "1"; exon_number "12"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding"; protein_id "ENSSHAP00000022389"; protein_version "1";
    # only keep between
    #@Devil txt=  'gene_id \"ENSSHAG00000011277\"; gene_version \"1\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\";'
    #@Devil txt <- gsub("gene_id.\"([^.]+)\"[;].*", "\\1", txt) #
    # now haveENSSHAG00000011277\"; gene_version \"1\"; gene_source \"ensembl\"; gene_biotype \"protein_coding
    # now we can split at the first \
    # strsplit(txt, '["]')[[1]][1] # ... in place
    # with our real data
    #@Devil txt <- gsub("gene_id.\"([^.]+)\"[;].*", "\\1", the.gtf.subset.TMP$attribute) #
    # loop
    #@Devil   try( remove(new.attribute.data)   )
    #@Devil   new.attribute.data <- c()
    #
    #@Devil   for (i in 1:(length(txt))) {
    #@Devil       the.data.we.want.for.this.sample <- txt[i]
    #@Devil       the.data.we.want.for.this.sample <- gsub("gene_id.\"([^.]+)\"[;].*", "\\1", the.data.we.want.for.this.sample) #
    #@Devil       the.data.we.want.for.this.sample <- strsplit(the.data.we.want.for.this.sample, '["]')[[1]][1]
    #@Devil       class(the.data.we.want.for.this.sample) # character
    #
    #@Devil       new.attribute.data[i] <- the.data.we.want.for.this.sample
    # cat(the.data.we.want.for.this.sample)
    #@Devil   }
    #  new.attribute.data
    
    
    
    
    
    
    
    
    
    # Antechinus flavipes
    #  [1] "transcript_id \"evm.model.ctg551_pilon_pilon.1\";"
    #
    # only keep between
    # txt=  'transcript_id \"evm.model.ctg551_pilon_pilon.1\";'
    # txt <- gsub('transcript_id'," ",txt)
    
    # now have ENSSHAG00000011277\"; gene_version \"1\"; gene_source \"ensembl\"; gene_biotype \"protein_coding
    # now we can split at the first \
    #    strsplit(txt, '["]')[[1]][1] # ... in place
    # with our real data
    
    
    txt <- gsub('transcript_id'," ", the.gtf.subset.TMP$attribute) #
    # loop
    try( remove(new.attribute.data)   )
    new.attribute.data <- c()
    for (i in 1:(length(txt))) {
        the.data.we.want.for.this.sample <- txt[i]
        the.data.we.want.for.this.sample <- gsub("transcript_id.\"([^.]+)\"[;].*", "\\1", the.data.we.want.for.this.sample) #
        the.data.we.want.for.this.sample <- strsplit(the.data.we.want.for.this.sample, '["]')[[1]][2]
        class(the.data.we.want.for.this.sample) # character
        #
        new.attribute.data[i] <- the.data.we.want.for.this.sample
        # cat(the.data.we.want.for.this.sample)
    }
    new.attribute.data
    # e.g. now have evm.model.ctg2_pilon_pilon.3, evm.model.ctg2_pilon_pilon.30, evm.model.ctg2_pilon_pilon.300 in dummy data
    
    # #####################################################################
    # to check that it is working, uncomment this line to only output a few ...
    #@ txt <- txt[1:10] # subset for now
    # #####################################################################
    
    
    
    
    
    
    
    # now have a list of all the genes
    # output a GFF file PER gene
    
    
    # next, keep only one gene ID
    one.of.each.gene <- new.attribute.data[!duplicated(new.attribute.data)]
    
    # so, in our 10 attribute example, we have two genes
    # ENSSHAG00000011277 ENSSHAG00000015543
    #@ i <- 2 # dummy
    
    library(data.table)
    
#    for (i in 1:(length(new.attribute.data))) {
#       desiredgene <- new.attribute.data[i]
        
        for (the.gene in 1:(length(one.of.each.gene))) {
            desiredgene <- one.of.each.gene[the.gene]
        
        
        
            
        #
        # temp.df <- the.gtf.subset.TMP[the.gtf.subset.TMP$attribute %like% desiredgene, ]
#        desired.gene 
#       test <- paste('"transcript_id \"',desiredgene,'\";"',sep="")
#        test <- paste('"transcript_id "',desiredgene,'";"',sep="")
        
        
#            temp.df <- the.gtf.subset.TMP[(grep(paste("evm.model.ctg2_pilon_pilon.3",'";', sep=""), the.gtf.subset.TMP$attribute)), ]
temp.df <- the.gtf.subset.TMP[(grep(paste(desiredgene,'";', sep=""), the.gtf.subset.TMP$attribute)), ]
            
        temp.file.name <- paste("./singleGTFs/",desiredgene,".gtf",sep="")
        #   class(temp.df) # df
        #    class(temp.file.name) # character
        write.table(temp.df, temp.file.name, sep = "\t",col.names = FALSE,quote = FALSE,row.names = FALSE,append = FALSE)
    }
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    
    
    
    
    
    
    
    
} # END OF BIG LOOP


# in the case of the gene ENSSHAG00000015543
# there two isoforms ... we will write each to a separate file later
#
# ATP11B-201    ENSSHAT00000018461.1    3843    1218aa
# ATP11B-202    ENSSHAT00000018462.1    3411    1137aa
# UniProt G3WSA7



# _THIS IS THE END

