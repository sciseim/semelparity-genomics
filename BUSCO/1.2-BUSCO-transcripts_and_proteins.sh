#!/usr/bin/env bash



# multiBUSCO

CPU=18

# STARTDIR

# need to edit config.ini # ... add BLAST etc. Cmd lines variables below will overrule some of these
# ERROR   No section [busco] found in /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/../config/config.ini. Please make sure both the file and this section exist, see userguide.


# BUSCO config
export AUGUSTUS_CONFIG_PATH=/bix/augustus-3.2.1/config/ # see https://busco.ezlab.org/v1/files/README.html
cd /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/pseudo/TasDev 
export BUSCO_CONFIG_FILE="${PWD}/config.ini"
head $BUSCO_CONFIG_FILE


# e.g. 
# n or X does not work... try with gaps
# Antechinus_argentus-TasDevPseudo-CDS.fa
#
# -m or --mode sets the assessment MODE: genome, proteins, transcriptome



# Gene set (proteins) assessment
cd /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/pseudo/TasDev/Protein
#
#
SAMPLENAME=Antechinus_argentus-TasDevPseudo-PROTEIN.fa # 
mkdir -p ./tmp 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# python scripts/run_BUSCO.py -i SEQUENCE_FILE -o OUTPUT_NAME -l LINEAGE -m prot
time python2.7 /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/run_BUSCO.py -i ${SAMPLENAME} --out ${SAMPLENAME}.protein_BUSCO --lineage_path /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/mammalia_odb9 --mode prot --cpu $CPU --tmp ${PWD}/tmp
# create graphic summary for BUSCO runs based on short summary files
cd run_${SAMPLENAME}.protein_BUSCO && python2.7 /bix/generate_plot.py -wd .
cd .. ;
#
SAMPLENAME=Antechinus_arktos-TasDevPseudo-PROTEIN.fa # 
mkdir -p ./tmp 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Gene set (proteins) assessment
# python scripts/run_BUSCO.py -i SEQUENCE_FILE -o OUTPUT_NAME -l LINEAGE -m prot
time python2.7 /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/run_BUSCO.py -i ${SAMPLENAME} --out ${SAMPLENAME}.protein_BUSCO --lineage_path /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/mammalia_odb9 --mode prot --cpu $CPU --tmp ${PWD}/tmp
# create graphic summary for BUSCO runs based on short summary files
cd run_${SAMPLENAME}.protein_BUSCO && python2.7 /bix/generate_plot.py -wd .
cd .. ;
#
SAMPLENAME=Murexia_melanurus-TasDevPseudo-PROTEIN.fa # 
mkdir -p ./tmp 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Gene set (proteins) assessment
# python scripts/run_BUSCO.py -i SEQUENCE_FILE -o OUTPUT_NAME -l LINEAGE -m prot
time python2.7 /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/run_BUSCO.py -i ${SAMPLENAME} --out ${SAMPLENAME}.protein_BUSCO --lineage_path /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/mammalia_odb9 --mode prot --cpu $CPU --tmp ${PWD}/tmp
# create graphic summary for BUSCO runs based on short summary files
cd run_${SAMPLENAME}.protein_BUSCO && python2.7 /bix/generate_plot.py -wd .
cd .. ;
#
SAMPLENAME=Murexia_sp-TasDevPseudo-PROTEIN.fa # 
mkdir -p ./tmp 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Gene set (proteins) assessment
# python scripts/run_BUSCO.py -i SEQUENCE_FILE -o OUTPUT_NAME -l LINEAGE -m prot
time python2.7 /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/run_BUSCO.py -i ${SAMPLENAME} --out ${SAMPLENAME}.protein_BUSCO --lineage_path /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/mammalia_odb9 --mode prot --cpu $CPU --tmp ${PWD}/tmp
# create graphic summary for BUSCO runs based on short summary files
cd run_${SAMPLENAME}.protein_BUSCO && python2.7 /bix/generate_plot.py -wd .
cd .. ;
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# REFERENCE
# cd /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/pseudo/TasDev/Protein ;
SAMPLENAME=TasDev # 
mkdir -p ./tmp 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Gene set (proteins) assessment
# python scripts/run_BUSCO.py -i SEQUENCE_FILE -o OUTPUT_NAME -l LINEAGE -m prot
time python2.7 /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/run_BUSCO.py -i /mnt/beakerdisk/ANTECHINUS/WGS/short-read/reference-based/TassieDevil/refgenome/output/AA.mfa --out ${SAMPLENAME}.protein_BUSCO --lineage_path /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/mammalia_odb9 --mode prot --cpu $CPU --tmp ${PWD}/tmp
# create graphic summary for BUSCO runs based on short summary files
cd run_${SAMPLENAME}.protein_BUSCO && python2.7 /bix/generate_plot.py -wd .
cd .. ;
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@









# 16.05.19
# does not appear to work on my CDS ... stick to the proteins. Way faster too

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Transcriptome assessment
# python scripts/run_BUSCO.py -i SEQUENCE_FILE -o OUTPUT_NAME -l LINEAGE -m prot
# Gene set (proteins) assessment
#@ cd /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/pseudo/TasDev/CDS
#
#
SAMPLENAME=Antechinus_argentus-TasDevPseudo-CDS.fa
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# python scripts/run_BUSCO.py -i SEQUENCE_FILE -o OUTPUT_NAME -l LINEAGE -m prot
time python2.7 /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/run_BUSCO.py -i ${SAMPLENAME} --out ${SAMPLENAME}.transcriptome_BUSCO --lineage_path /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/mammalia_odb9 --mode transcriptome --cpu $CPU --tmp ${PWD}/tmp
# create graphic summary for BUSCO runs based on short summary files
cd run_${SAMPLENAME}.transcriptome_BUSCO && python2.7 /bix/generate_plot.py -wd .
cd .. ;
#
SAMPLENAME=Antechinus_arktos-TasDevPseudo-CDS.fa
.fa # 
mkdir -p ./tmp 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Gene set (proteins) assessment
# python scripts/run_BUSCO.py -i SEQUENCE_FILE -o OUTPUT_NAME -l LINEAGE -m prot
time python2.7 /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/run_BUSCO.py -i ${SAMPLENAME} --out ${SAMPLENAME}.transcriptome_BUSCO --lineage_path /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/mammalia_odb9 --mode transcriptome
 --cpu $CPU --tmp ${PWD}/tmp
# create graphic summary for BUSCO runs based on short summary files
cd run_${SAMPLENAME}.transcriptome_BUSCO && python2.7 /bix/generate_plot.py -wd .
cd .. ;
#
SAMPLENAME=Murexia_melanurus-TasDevPseudo-CDS.fa
.fa # 
mkdir -p ./tmp 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Gene set (proteins) assessment
# python scripts/run_BUSCO.py -i SEQUENCE_FILE -o OUTPUT_NAME -l LINEAGE -m prot
time python2.7 /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/run_BUSCO.py -i ${SAMPLENAME} --out ${SAMPLENAME}.transcriptome_BUSCO --lineage_path /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/mammalia_odb9 --mode transcriptome
 --cpu $CPU --tmp ${PWD}/tmp
# create graphic summary for BUSCO runs based on short summary files
cd run_${SAMPLENAME}.transcriptome_BUSCO && python2.7 /bix/generate_plot.py -wd .
cd .. ;
#
SAMPLENAME=Murexia_sp-TasDevPseudo-CDS.fa
.fa # 
mkdir -p ./tmp 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Gene set (proteins) assessment
# python scripts/run_BUSCO.py -i SEQUENCE_FILE -o OUTPUT_NAME -l LINEAGE -m prot
time python2.7 /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/run_BUSCO.py -i ${SAMPLENAME} --out ${SAMPLENAME}.transcriptome_BUSCO --lineage_path /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/mammalia_odb9 --mode transcriptome
 --cpu $CPU --tmp ${PWD}/tmp
# create graphic summary for BUSCO runs based on short summary files
cd run_${SAMPLENAME}.transcriptome_BUSCO && python2.7 /bix/generate_plot.py -wd .
cd .. ;
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# REFERENCE
# cd /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/pseudo/TasDev/Protein ;
SAMPLENAME=TasDev # 
mkdir -p ./tmp 
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Gene set (proteins) assessment
# python scripts/run_BUSCO.py -i SEQUENCE_FILE -o OUTPUT_NAME -l LINEAGE -m prot
time python2.7 /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/run_BUSCO.py -i /mnt/beakerdisk/ANTECHINUS/WGS/short-read/reference-based/TassieDevil/refgenome/output/AA.mfa --out ${SAMPLENAME}.transcriptome_BUSCO --lineage_path /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/mammalia_odb9 --mode transcriptome --cpu $CPU --tmp ${PWD}/tmp
# create graphic summary for BUSCO runs based on short summary files
cd run_${SAMPLENAME}.transcriptome_BUSCO && python2.7 /bix/generate_plot.py -wd .
cd .. ;
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



# _THIS IS THE END