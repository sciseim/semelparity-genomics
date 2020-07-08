#!/usr/bin/env bash

# from excellent GitHub repo https://github.com/jwillou/wgs_compare/
cd /mnt/beakerdisk/ANTECHINUS/PSMC/mutation_rate/TasDevPseudo-vs-TasDev/wgs_compare

# create a BLASTdb of REFERENCE
makeblastdb -in /mnt/beakerdisk/ANTECHINUS/PSMC/mutation_rate/TasDevPseudo-vs-TasDev/Sarcophilus_harrisii.DEVIL7.0.dna_rm.toplevel.fa -dbtype nucl
touch todivR.txt


# for each pseudogenome
SAMPLENAME=AA100A
mkdir ${SAMPLENAME}_vs_TasDev
echo ${SAMPLENAME}_vs_TasDev >> todivR.txt # add sample info to file for later parsing
#
perl ./fastasplit_chunks1000.pl /mnt/beakerdisk/ANTECHINUS/PSMC/mutation_rate/TasDevPseudo-vs-TasDev/${SAMPLENAME}.unique_d10_Q30_genome.fasta >${SAMPLENAME}_1000.fa 
blastn -query ${SAMPLENAME}_1000.fa  -db /mnt/beakerdisk/ANTECHINUS/PSMC/mutation_rate/TasDevPseudo-vs-TasDev/Sarcophilus_harrisii.DEVIL7.0.dna_rm.toplevel.fa -out ./${SAMPLENAME}_vs_TasDev/outblast21.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident mismatch gapopen gaps qseq sseq' -perc_identity 3 -qcov_hsp_perc 3 -culling_limit 1 -max_target_seqs 1 -evalue 10 -num_threads 15
# @@@@@@@@@@


# for each pseudogenome
SAMPLENAME=BD-17-5A
mkdir ${SAMPLENAME}_vs_TasDev
echo ${SAMPLENAME}_vs_TasDev >> todivR.txt # add sample info to file for later parsing
#
perl ./fastasplit_chunks1000.pl /mnt/beakerdisk/ANTECHINUS/PSMC/mutation_rate/TasDevPseudo-vs-TasDev/${SAMPLENAME}.unique_d10_Q30_genome.fasta >${SAMPLENAME}_1000.fa 
blastn -query ${SAMPLENAME}_1000.fa  -db /mnt/beakerdisk/ANTECHINUS/PSMC/mutation_rate/TasDevPseudo-vs-TasDev/Sarcophilus_harrisii.DEVIL7.0.dna_rm.toplevel.fa -out ./${SAMPLENAME}_vs_TasDev/outblast21.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident mismatch gapopen gaps qseq sseq' -perc_identity 3 -qcov_hsp_perc 3 -culling_limit 1 -max_target_seqs 1 -evalue 10 -num_threads 15
# @@@@@@@@@@
# for each pseudogenome
SAMPLENAME=46020A
mkdir ${SAMPLENAME}_vs_TasDev
echo ${SAMPLENAME}_vs_TasDev >> todivR.txt # add sample info to file for later parsing
#
perl ./fastasplit_chunks1000.pl /mnt/beakerdisk/ANTECHINUS/PSMC/mutation_rate/TasDevPseudo-vs-TasDev/${SAMPLENAME}.unique_d10_Q30_genome.fasta >${SAMPLENAME}_1000.fa 
blastn -query ${SAMPLENAME}_1000.fa  -db /mnt/beakerdisk/ANTECHINUS/PSMC/mutation_rate/TasDevPseudo-vs-TasDev/Sarcophilus_harrisii.DEVIL7.0.dna_rm.toplevel.fa -out ./${SAMPLENAME}_vs_TasDev/outblast21.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident mismatch gapopen gaps qseq sseq' -perc_identity 3 -qcov_hsp_perc 3 -culling_limit 1 -max_target_seqs 1 -evalue 10 -num_threads 15
# @@@@@@@@@@
# for each pseudogenome
SAMPLENAME=42595A
mkdir ${SAMPLENAME}_vs_TasDev
echo ${SAMPLENAME}_vs_TasDev >> todivR.txt # add sample info to file for later parsing
#
perl ./fastasplit_chunks1000.pl /mnt/beakerdisk/ANTECHINUS/PSMC/mutation_rate/TasDevPseudo-vs-TasDev/${SAMPLENAME}.unique_d10_Q30_genome.fasta >${SAMPLENAME}_1000.fa 
blastn -query ${SAMPLENAME}_1000.fa  -db /mnt/beakerdisk/ANTECHINUS/PSMC/mutation_rate/TasDevPseudo-vs-TasDev/Sarcophilus_harrisii.DEVIL7.0.dna_rm.toplevel.fa -out ./${SAMPLENAME}_vs_TasDev/outblast21.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident mismatch gapopen gaps qseq sseq' -perc_identity 3 -qcov_hsp_perc 3 -culling_limit 1 -max_target_seqs 1 -evalue 10 -num_threads 15
# @@@@@@@@@@



# for each pseudogenome
SAMPLENAME=Koala
mkdir ${SAMPLENAME}_vs_TasDev
echo ${SAMPLENAME}_vs_TasDev >> todivR.txt # add sample info to file for later parsing
#
# perl ./fastasplit_chunks1000.pl /mnt/beakerdisk/ANTECHINUS/ROPUS/genomes_ncbi/phascolarctos_cinereus/phascolarctos_cinereus_ncbi.fa >${SAMPLENAME}_1000.fa 
# @@@@@@@@@@
# perl ./fastasplit_chunks1000.pl /mnt/beakerdisk/ANTECHINUS/ROPUS/genomes_ensembl/phascolarctos_cinereus/phascolarctos_cinereus_ensembl.fa >${SAMPLENAME}_1000.fa 
# for some reason it does NOT for koal (1900 or so sequences... this splits into 1000 bp chunks)
python ./splitter.py /mnt/beakerdisk/ANTECHINUS/ROPUS/genomes_ncbi/phascolarctos_cinereus/phascolarctos_cinereus_ncbi.fa 1000 >>${SAMPLENAME}_1000.fa 
blastn -query ${SAMPLENAME}_1000.fa  -db /mnt/beakerdisk/ANTECHINUS/PSMC/divergence_rate/TasDevPseudo-vs-TasDev/BLASTdb/Sarcophilus_harrisii.DEVIL7.0.dna_rm.toplevel.fa -out ./${SAMPLENAME}_vs_TasDev/outblast21.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident mismatch gapopen gaps qseq sseq' -perc_identity 3 -qcov_hsp_perc 3 -culling_limit 1 -max_target_seqs 1 -evalue 10 -num_threads 15

# AntFla
SAMPLENAME=AntFla
mkdir ${SAMPLENAME}_vs_TasDev
echo ${SAMPLENAME}_vs_TasDev >> todivR.txt # add sample info to file for later parsing
#
python ./splitter.py /mnt/beakerdisk/ANTECHINUS/WGS/BUSCO/AdamAnt/Yellowfoot.map.fa 1000 >>${SAMPLENAME}_1000.fa 
blastn -query ${SAMPLENAME}_1000.fa  -db /mnt/beakerdisk/ANTECHINUS/PSMC/divergence_rate/TasDevPseudo-vs-TasDev/BLASTdb/Sarcophilus_harrisii.DEVIL7.0.dna_rm.toplevel.fa -out ./${SAMPLENAME}_vs_TasDev/outblast21.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident mismatch gapopen gaps qseq sseq' -perc_identity 3 -qcov_hsp_perc 3 -culling_limit 1 -max_target_seqs 1 -evalue 10 -num_threads 15

# AntFlaPseudo
SAMPLENAME=AntFlaPseudo
mkdir ${SAMPLENAME}_vs_TasDev
echo ${SAMPLENAME}_vs_TasDev >> todivR.txt # add sample info to file for later parsing
#
python ./splitter.py /mnt/beakerdisk/ANTECHINUS/PSMC/pseudogenomes-no-chrX/flavipes/flavipes.unique_d10_Q30_genome.fasta 1000 >>${SAMPLENAME}_1000.fa 
blastn -query ${SAMPLENAME}_1000.fa  -db /mnt/beakerdisk/ANTECHINUS/PSMC/divergence_rate/TasDevPseudo-vs-TasDev/BLASTdb/Sarcophilus_harrisii.DEVIL7.0.dna_rm.toplevel.fa -out ./${SAMPLENAME}_vs_TasDev/outblast21.txt -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore nident mismatch gapopen gaps qseq sseq' -perc_identity 3 -qcov_hsp_perc 3 -culling_limit 1 -max_target_seqs 1 -evalue 10 -num_threads 15
