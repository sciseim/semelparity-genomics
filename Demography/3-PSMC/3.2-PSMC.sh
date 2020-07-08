#!/usr/bin/env bash

# ##################################################################

# run this after pseudo-genome-creator.sh


SAMPLENAME=46020A
MAXREADDEPTH=55.6 # twice the average genome coverage
PSEUDOGENOME=/mnt/beakerdisk/ANTECHINUS/PSMC/pseudogenomes-no-chrX/${SAMPLENAME}/${SAMPLENAME}.unique_d10_Q30_genome.fasta
PSEUDOBAM=/mnt/beakerdisk/ANTECHINUS/PSMC/pseudogenomes-no-chrX/${SAMPLENAME}/${SAMPLENAME}.aln.sorted.unique.bam

# SAMPLENAME=42595A
OLDSAMTOOLSPATH=/mnt/beakerdisk/ANTECHINUS/WGS/short-read/reference-based/TassieDevil/RefGenBin/samtools-0.1.19
BCFTOOLSPATH=${OLDSAMTOOLSPATH}/bcftools/


# from https://informatics.fas.harvard.edu/psmc-journal-club-walkthrough.html
# samtools:
# -Q and -q in mpileup determine the cutoffs for baseQ and mapQ, respectively
# -v tells mpileup to produce vcf output, and -u says that should be uncompressed
# -f is the reference fasta used (needs to be indexed)
# -r is the region to call the mpileup for (in this case, a particular chromosome based on the array task id)
# P964.bam is the bam file to use
# bcftools:
# call -c calls a consensus sequence from the mpileup using the original calling method
# vcfutils.pl:
# -d 5 and -d 34 determine the minimum and maximum coverage to allow for vcf2fq, anything outside that range is filtered
# -Q 30 sets the root mean squared mapping quality minimum to 30

# 080419 -- run newer samtools as per HMS tutorial
# THEBIGCOMMAND=$($OLDSAMTOOLSPATH/samtools mpileup -Q 30 -q 30 -u -v -f $REFERENCEGENOMEPATH ${SAMPLENAME}.aln.sorted.unique.bam | ${BCFTOOLSPATH}/bcftools call -c | ${BCFTOOLSPATH}/vcfutils.pl vcf2fq -d 5 -D 34 -Q 30 > ${SAMPLENAME}_for_PSMC.fq)
# time $THEBIGCOMMAND


# We'll use Odyssey's module system to control software versions. For reference, we use samtools version 1.2 and bcftools version 1.2. Note that these are newer versions than used in the paper.
# v1.7
# THEBIGCOMMAND=$(samtools mpileup -Q 30 -q 30 -u -v -f $REFERENCEGENOMEPATH ${SAMPLENAME}.aln.sorted.unique.bam | bcftools call -c | ${BCFTOOLSPATH}/vcfutils.pl vcf2fq -d 5 -D 34 -Q 30 > ${SAMPLENAME}_for_PSMC.fq)
# time $THEBIGCOMMAND


# ###################################################################
# map FASTQ files to a diploid genome and obtain heterozygous SNPs
# ###################################################################

# actually, I want to run it against the pseudo-genome!
# PSEUDOGENOME=${SAMPLENAME}.unique_d10_Q30_genome.fasta
samtools mpileup -Q 30 -q 30 -u -v -f $PSEUDOGENOME $PSEUDOBAM | bcftools call -c | ${BCFTOOLSPATH}/vcfutils.pl vcf2fq -d 5 -D ${MAXREADDEPTH} -Q 30 > ${SAMPLENAME}_for_PSMC.fq


# samtools mpileup -Q 30 -q 30 -u -v -f $PSEUDOGENOME ${SAMPLENAME}.aln.sorted.unique.bam | bcftools call -c | ${BCFTOOLSPATH}/vcfutils.pl vcf2fq -d 5 -D 34 -Q 30 > ${SAMPLENAME}_for_PSMC.fq


# THEBIGCOMMAND=$(samtools mpileup -Q 30 -q 30 -u -v -f $PSEUDOGENOME ${SAMPLENAME}.aln.sorted.unique.bam | bcftools call -c | ${BCFTOOLSPATH}/vcfutils.pl vcf2fq -d 5 -D 34 -Q 30 > ${SAMPLENAME}_for_PSMC.fq)
# time $THEBIGCOMMAND




# ###################################################################
# PSMC
# ###################################################################

# Infer the history of population sizes
# https://github.com/lh3/psmc
# The pairwise sequentially Markovian coalescent (PSMC) method uses the genome sequence of a single individual to estimate demographic history covering a time span of thousands of generations
PSMCPATH=/mnt/beakerdisk/ANTECHINUS/WGS/short-read/reference-based/TassieDevil/RefGenBin/psmc-master


# perform bootstrapping
# To perform bootstrapping, one has to run splitfa first to split long chromosome
# sequences to shorter segments. When the `-b' option is applied, psmc will then
# randomly sample with replacement from these segments. As an example, the
# following command lines perform 100 rounds of bootstrapping:
${PSMCPATH}/utils/fq2psmcfa -q20 ${SAMPLENAME}_for_PSMC.fq > ${SAMPLENAME}_for_PSMC.psmcfa

# Split for bootstrapping
${PSMCPATH}/utils/splitfa ${SAMPLENAME}_for_PSMC.psmcfa > ${SAMPLENAME}_split.psmcfa

# Run PSMC on the genome
${PSMCPATH}/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o ${SAMPLENAME}.psmc ${SAMPLENAME}_for_PSMC.psmcfa


# later bootstrap! remember -- must be run in parallel!!


