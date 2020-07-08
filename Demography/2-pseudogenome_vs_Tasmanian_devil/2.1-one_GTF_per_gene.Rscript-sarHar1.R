# to run: from Rscript ./one_GTF_per_gene.Rscript $REFERENCEGTFPATH
args <- commandArgs()
print(args)

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
# it handles doublequote automatically e.g. #!

# colnames(the.gtf)
# the.gtf[1:1,7:9]
# head(the.gtf)
# class(the.gtf)


# dummy until I am 100% sure it is working
# the.gtf <- the.gtf[1:10000,] # Testing: 347 genes
# tail(the.gtf)
# ENSSHAG00000016243 ... should be at least three genes now... is 5




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
length(unique.chr) # 15 scaffolds

# i <- 3 # dummy 
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
  # [1000] gene_id "ENSSHAG00000018951"; gene_version "1"; transcript_id "ENSSHAT00000022567"; transcript_version "1"; exon_number "12"; gene_source "ensembl"; gene_biotype "protein_coding"; transcript_source "ensembl"; transcript_biotype "protein_coding"; protein_id "ENSSHAP00000022389"; protein_version "1";  
  
  # only keep between
  txt=  'gene_id \"ENSSHAG00000011277\"; gene_version \"1\"; gene_source \"ensembl\"; gene_biotype \"protein_coding\";'
  txt <- gsub("gene_id.\"([^.]+)\"[;].*", "\\1", txt) #
  # now haveENSSHAG00000011277\"; gene_version \"1\"; gene_source \"ensembl\"; gene_biotype \"protein_coding
  # now we can split at the first \
  strsplit(txt, '["]')[[1]][1] # ... in place
  
  # with our real data
  txt <- gsub("gene_id.\"([^.]+)\"[;].*", "\\1", the.gtf.subset.TMP$attribute) #
  # loop
  
  
  # #####################################################################
  # to check that it is working, uncomment this line to only output a few ...
  #@ txt <- txt[1:10] # subset for now
  # #####################################################################
  
  
  
  
  
  try( remove(new.attribute.data)   )
  new.attribute.data <- c()
  #
  for (i in 1:(length(txt))) {
    the.data.we.want.for.this.sample <- txt[i]
    the.data.we.want.for.this.sample <- gsub("gene_id.\"([^.]+)\"[;].*", "\\1", the.data.we.want.for.this.sample) #
    the.data.we.want.for.this.sample <- strsplit(the.data.we.want.for.this.sample, '["]')[[1]][1]
    class(the.data.we.want.for.this.sample) # character
    #
    new.attribute.data[i] <- the.data.we.want.for.this.sample
    # cat(the.data.we.want.for.this.sample)
  }
#  new.attribute.data
  
  # now have a list of all the genes 
  # output a GFF file PER gene
  
  
  # next, keep only one gene ID
  one.of.each.gene <- new.attribute.data[!duplicated(new.attribute.data)]
  
  # so, in our 10 attribute example, we have two genes
  # ENSSHAG00000011277 ENSSHAG00000015543
  #@ i <- 2 # dummy
  
  library(data.table)
  
  for (i in 1:(length(new.attribute.data))) {
    desiredgene <- new.attribute.data[i]
    
 
    #  
    temp.df <- the.gtf.subset.TMP[the.gtf.subset.TMP$attribute %like% desiredgene, ]
    temp.file.name <- paste("./singleGTFs/",desiredgene,".gtf",sep="")
 #   class(temp.df) # df
 #    class(temp.file.name) # character
    write.table(temp.df, temp.file.name, sep = "\t",col.names = FALSE,quote = FALSE,row.names = FALSE,append = FALSE)   
  }  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
  
  
  
  
  
  
  
} # END OF BIG LOOP


# in the case of the gene ENSSHAG00000015543
# there are two isoforms ... we will write each to a separate file later
# 
# ATP11B-201	ENSSHAT00000018461.1	3843	1218aa	
# ATP11B-202	ENSSHAT00000018462.1	3411	1137aa	
# UniProt G3WSA7	



# _THIS IS THE END