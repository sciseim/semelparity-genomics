# #!/usr/bin/env bash

# pseudo-genome-creator.sh

# #############################################################
# Pseudo-genome pipeline
#
# by Inge Seim 
# ... last updated 06.04.2019
# #############################################################


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# NOTE: need to repeat-mask the Antechinus flavipes genome before mapping!
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# logic behind this approach
# used by
# Genome of the Tasmanian tiger provides insights into the evolution and demography of an extinct marsupial carnivore   :  https://www.nature.com/articles/s41559-017-0417-y
# Whole-genome sequencing of the blue whale and other rorquals finds signatures for introgressive gene flow : http://advances.sciencemag.org/content/4/4/eaap9873
# Complete genomes reveal signatures of demographic and genetic declines in the woolly mammoth: https://www.ncbi.nlm.nih.gov/pubmed/25913407

# Gaps/N are introduced whereever other species do not have reads matching A. flavipes, so can use A. flavipes genome annotation (gff) to pull out the other species' genes/coding sequences!
# [gene A].....NNNNN (masked REPEATS)......[gene Y]
#   ATGC                                      ATGC    Antechinus flavipes reference genome
#   ATGA                                      CNNC    Other, related species
#
# ... but first we must generate a gapfilled reference-based genome of our species of interest. E.g. Murexia sp. vs Tasmanian devil







# Finished
# 42595A

# to-do
# BD-17-5A
# /mnt/beakerdisk/ANTECHINUS/WGS/short-read/reference-based/TassieDevil/BD-17-5A/BD-17-5A_uncontaminated_noMT.1.fastq
# /mnt/beakerdisk/ANTECHINUS/WGS/short-read/reference-based/TassieDevil/BD-17-5A/BD-17-5A_uncontaminated_noMT.2.fastq
# 46020A
# /mnt/beakerdisk/ANTECHINUS/WGS/short-read/reference-based/TassieDevil/46020A/46020A_uncontaminated_noMT.1.fastq
# /mnt/beakerdisk/ANTECHINUS/WGS/short-read/reference-based/TassieDevil/46020A/46020A_uncontaminated_noMT.2.fastq

# AA100A
# /mnt/beakerdisk/ANTECHINUS/WGS/short-read/reference-based/TassieDevil/AA100A/AA100A_uncontaminated_noMT.1.fastq
# /mnt/beakerdisk/ANTECHINUS/WGS/short-read/reference-based/TassieDevil/AA100A/AA100A_uncontaminated_noMT.2.fastq




# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
# custom settings here

# edit this for each sample
CPU=18  # set to the number of CPUs at hand
SAMPLENAME=46020A
READ1=/mnt/beakerdisk/ANTECHINUS/WGS/short-read/reference-based/TassieDevil/46020A/46020A_uncontaminated_noMT.1.fastq
READ2=/mnt/beakerdisk/ANTECHINUS/WGS/short-read/reference-based/TassieDevil/46020A/46020A_uncontaminated_noMT.2.fastq
# head $READ1
# head $READ2

# IMPORTANT IMPORTANT IMPORTANT IMPORTANT
# the reference genome FASTA *must* be repeat-masked!
# IMPORTANT IMPORTANT IMPORTANT IMPORTANT
#
# complexity regions are detected with the RepeatMasker tool and masked
# by replacing repeats with 'N's.
#
# # note: ensembl has masked the Tassie devil genome (same as NCBI GCA_000189315.1) for us!  'dna_rm' - masked genomic DNA.  
REFERENCEGENOMEPATH=/mnt/beakerdisk/ANTECHINUS/PSMC/pseudogenomes-no-chrX/putative_TasDev_chrX-linked_scaffolds/TasDev_no_chrX.fa
REFERENCEGTFPATH=/mnt/beakerdisk/ANTECHINUS/WGS/short-read/reference-based/TassieDevil/refgenome/TasDev/Sarcophilus_harrisii.DEVIL7.0.95.gtf



# sometimes a GTF entry is 'difficult', breaking the parsing, so split your reference GTF file into individual files
# use 
# one_GTF_per_gene.Rscript.R
# to run: from Rscript ./one_GTF_per_gene.Rscript $REFERENCEGTFPATH
# ... OBVIOUSLY ONLY NEED TO RUN IT ONCE ... 
# it will put the reference genome's GTF files in ./gtf
SINGLEGTFSPATH=/mnt/beakerdisk/ANTECHINUS/WGS/short-read/reference-based/TassieDevil/refgenome/TasDev/singleGTFs



# remember to INDEX the genome for BWA too
# bwa index nameofrefgenome.fa # in the same folder as the reference genome






# /\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\

# -----------------------------------------------
# the following tools are required
# 040419: NEWER versions of these tool do not work! Stick to the versions below so we can use parameters from various pseudo-genome papers
# bwa mem v0.7.12
# samtools v0.1.19
# bcftools v0.1.19 # IN samtools!
# vcfutils # a Perl script in bcftools
# USE THESE OLDER VERSIONS ! store them in ./RefGenBin/
#
# tar xjf bwa-0.7.12.tar.bz2
# cd bwa-0.7.12 && make
# to run this older version
OLDBWA=/mnt/beakerdisk/ANTECHINUS/WGS/short-read/reference-based/TassieDevil/RefGenBin/bwa-0.7.12/bwa
echo $OLDBWA
#
# tar xjf samtools-0.1.19.tar.bz2
# cd ../samtools-0.1.19 && make
# cd ./bcftools && make
OLDSAMTOOLSPATH=/mnt/beakerdisk/ANTECHINUS/WGS/short-read/reference-based/TassieDevil/RefGenBin/samtools-0.1.19
BCFTOOLSPATH=${OLDSAMTOOLSPATH}/bcftools/
# -----------------------------------------------


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Map fastq reads to REFERENCE genome (e.g. Tasmanian devil or Antechinus flavipes)
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# BWA-MEM is a new alignment algorithm for aligning sequence reads or long query sequences against a large reference genome such as human. It automatically chooses between local and end-to-end alignments, supports paired-end reads and performs chimeric alignment. The algorithm is robust to sequencing errors and applicable to a wide range of sequence lengths from 70bp to a few megabases. For mapping 100bp sequences, BWA-MEM shows better performance than several state-of-art read aligners to date.
# https://github.com/jts/sga/issues/121

# bwa mem -t $CPU $REFERENCEGENOMEPATH $READ1 $READ2 | samtools view -F2304 -b -o reads.bam -
#       -B INT        penalty for a mismatch [4]
#       -O INT[,INT]  gap open penalties for deletions and insertions [6,6]
# map reads to repeat-masked REFERENCE genome
time $OLDBWA mem -B 3 -O 5,5 -t $CPU $REFERENCEGENOMEPATH $READ1 $READ2 > ${SAMPLENAME}.aln.sam
# run time: 411m ~7h
# 273 GB SAM file

# Converting SAM to BAM
time $OLDSAMTOOLSPATH/samtools view -bS ${SAMPLENAME}.aln.sam > ${SAMPLENAME}.aln.bam ;
# run time: 185m=3h
# 100.7 GB BAM file

# can probably delete the samfile to save space = at least ~350 GB
rm ${SAMPLENAME}.aln.sam

# Sorting and Indexing
# As our goal is to call genomic variants, and this requires that we “pile-up” all matching reads within a specific genomic location, we sort by location:
time $OLDSAMTOOLSPATH/samtools sort ${SAMPLENAME}.aln.bam ${SAMPLENAME}.aln.sorted
# run time: 233m ~4h
# 77.2 GB .aln.sorted.BAM file

# can probably delete the unsorted BAM file to save space = at least ~100 GB
rm ${SAMPLENAME}.aln.bam

# Once you have sorted your BAM file, you can then index it. This enables tools, including SAMtools itself, and other genomic viewers to perform efficient random access on the BAM file, resulting in greatly improved performance. To do so, run:
time $OLDSAMTOOLSPATH/samtools index ${SAMPLENAME}.aln.sorted.bam
# run time: 12m
# 9.7 MB index file


# remove PCR duplicates
# Duplicate filtering may improve predictive accuracy relative to no filtering
# https://www.ncbi.nlm.nih.gov/pubmed/27454357
# ...Our results suggest that PCR duplicate removal has minimal effect on the accuracy of subsequent variant calls.
# ...Picard required both more memory and execution time than SAMTools
# ... also a must for coverage stats, so let us just do it ...

# Potential PCR duplicates within each set of reads having a common molecular index were removed with samtools rmdup.
# the normal % of duplicates appears to be about 5%
time $OLDSAMTOOLSPATH/samtools rmdup ${SAMPLENAME}.aln.sorted.bam ${SAMPLENAME}.aln.sorted.unique.bam


# can probably delete the BAM file to save space = at least ~70 GB
rm ${SAMPLENAME}.aln.sorted.bam
rm ${SAMPLENAME}.aln.sorted.bai






# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# pile-up” all matching reads within a specific genomic location
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# create reference-based assembly
# mpileup        multi-way pileup
# 
# -B, --no-BAQ            disable BAQ (per-Base Alignment Quality)
# -U FILE Write alignments that are not selected by the various filter options to FILE. When this option is used, all alignments (or all alignments intersecting the regions specified) are written to either the output file or this file, but never both.
# -I, --id STR
# Include only listed read group or sample name []
# Usage: samtools mpileup [options] in1.bam [in2.bam [...]]

# parameters −I −B −u, without inclusion of the reference genome fasta
# KEY: do NOT print-out the bases of the REFERENCE species -- ONLY your short-insert WGS species (e.g. Murexia sp.)!!!


# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# create a VCF file using samtools
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
time $OLDSAMTOOLSPATH/samtools mpileup -I -B -u ${SAMPLENAME}.aln.sorted.unique.bam | ${BCFTOOLSPATH}/bcftools view -I -cg - >${SAMPLENAME}.unique.vcf ;


# ####################################################
# convert to FASTQ
# ####################################################
# what they mean ... default indicated
# -Q INT     minimum RMS mapping quality for SNPs [10]
# -d INT      minimum read depth [2]
# -d INT  maximum read depth [10000000]



# filtering
# woolly mammoth MS: 

# we can set various minimum read coverage, but would be wise to know what it is in our sample, rather than relying on published short-insert pseudo-genome MS
# average coverage?
# samtools depth bamfile | awk '{sum+=$3} END { print "Average = ",sum/NR}'
AVERAGECOVERAGE=$($OLDSAMTOOLSPATH/samtools depth ${SAMPLENAME}.aln.sorted.unique.bam | awk '{sum+=$3} END { print "Average = ",sum/NR}')
echo $AVERAGECOVERAGE
echo $AVERAGECOVERAGE > ${SAMPLENAME}.average_coverage.txt
# Average = 27.5262 # *good* -- expected about 30X
# 


# ####################################################
# convert FASTQ to FASTA
# ####################################################
# git clone https://github.com/lh3/misc.git
# cd misc/seq/seqtk/
# gcc -O2 seqtk.c -o seqtk -lz
# seqtk fq2fa in.fastq > out.fasta

# Convert FASTQ to FASTA and set bases of quality lower than 30 to N




# may stick with d=10 and Q=30 ... mammoth MS in Current Biology
time perl ${BCFTOOLSPATH}/vcfutils.pl vcf2fq -d 10 ${SAMPLENAME}.unique.vcf >${SAMPLENAME}.unique_d10_genome.fastq


# can probably delete the BAM file to save space = at least ~160 GB
# takes longer to retrieve from Synology server than re-running the above step!
rm ${SAMPLENAME}.unique.vcf




time seqtk seq -aQ64 -q30 -n N ${SAMPLENAME}.unique_d10_genome.fastq >${SAMPLENAME}.unique_d10_Q30_genome.fasta  # Q64 here since BGISEQ-500

# ####################################################
# obtain pseudo-genome assembly statistics
# ####################################################
# of course, silly since we have mapped to another genome, but useful nevertheless to e.g. assess how many gaps, etc

# sudo apt install assembly-stats
# https://github.com/sanger-pathogens/assembly-stats

# input species
time assembly-stats ${SAMPLENAME}.unique_d10_Q30_genome.fasta > ${SAMPLENAME}_assembly_stats.txt

# the reference
# note: Tasmanian devil genome has quite a few unplaced scaffolds, etc.
time assembly-stats $REFERENCEGENOMEPATH > REFERENCE_assembly_stats.txt









# _THIS IS THE END

