#!/usr/bin/env bash

CPU=10
# STARTDIR
STARTDIR=/mnt/beakerdisk/ANTECHINUS/WGS/BUSCO/AdamAnt/HiC_for_MS/
# note, please set CPU=1 if you have a LARGE number of contigs/scaffolds in your genome assembly (BUSCO bug)

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# for each
cd $STARTDIR
# 42595A.unique_d10_Q30_genome.fasta
SAMPLENAME=AdamAnt
cd $SAMPLENAME
# need to edit config.ini # ... add BLAST etc. Cmd lines variables below will overrule some of these
# ERROR   No section [busco] found in /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/../config/config.ini. Please make sure both the file and this section exist, see userguide.
#
mkdir -p ./tmp 
export AUGUSTUS_CONFIG_PATH=/bix/augustus-3.2.1/config/ # see https://busco.ezlab.org/v1/files/README.html
export BUSCO_CONFIG_FILE="${PWD}/config.ini"
head $BUSCO_CONFIG_FILE
#

# sga
# python2.7 /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/run_BUSCO.py -i ${SAMPLENAME}_gapfilled.fa --out $SAMPLENAME --lineage_path /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/mammalia_odb9 --mode genome --cpu $CPU --tmp ${PWD}/tmp
python2.7 /mnt/beakerdisk/ANTECHINUS/WGS/BUSCO/busco-master/scripts/run_BUSCO.py -i ${SAMPLENAME}.fa --out $SAMPLENAME --lineage_path /mnt/beakerdisk/ANTECHINUS/WGS/BUSCO/mammalia_odb9 --mode genome --cpu $CPU --tmp ${PWD}/tmp
# Human genome (3.1 Gbp), assessed with 4’104 mammalian BUSCOs: 6 days 15 hours
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# for each
cd $STARTDIR
# 42595A.unique_d10_Q30_genome.fasta
SAMPLENAME=AA100A
cd $SAMPLENAME
# need to edit config.ini # ... add BLAST etc. Cmd lines variables below will overrule some of these
# ERROR   No section [busco] found in /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/../config/config.ini. Please make sure both the file and this section exist, see userguide.
#
mkdir -p ./tmp 
export AUGUSTUS_CONFIG_PATH=/bix/augustus-3.2.1/config/ # see https://busco.ezlab.org/v1/files/README.html
export BUSCO_CONFIG_FILE="${PWD}/config.ini"
head $BUSCO_CONFIG_FILE
#

# sga
# python2.7 /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/run_BUSCO.py -i ${SAMPLENAME}_gapfilled.fa --out $SAMPLENAME --lineage_path /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/mammalia_odb9 --mode genome --cpu $CPU --tmp ${PWD}/tmp
python2.7 /mnt/beakerdisk/ANTECHINUS/WGS/BUSCO/busco-master/scripts/run_BUSCO.py -i ${SAMPLENAME}.fa --out $SAMPLENAME --lineage_path /mnt/beakerdisk/ANTECHINUS/WGS/BUSCO/mammalia_odb9 --mode genome --cpu $CPU --tmp ${PWD}/tmp
# Human genome (3.1 Gbp), assessed with 4’104 mammalian BUSCOs: 6 days 15 hours
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# for each
cd $STARTDIR
# 42595A.unique_d10_Q30_genome.fasta
SAMPLENAME=BD-17-5A
cd $SAMPLENAME
# need to edit config.ini # ... add BLAST etc. Cmd lines variables below will overrule some of these
# ERROR   No section [busco] found in /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/../config/config.ini. Please make sure both the file and this section exist, see userguide.
#
mkdir -p ./tmp 
export AUGUSTUS_CONFIG_PATH=/bix/augustus-3.2.1/config/ # see https://busco.ezlab.org/v1/files/README.html
export BUSCO_CONFIG_FILE="${PWD}/config.ini"
head $BUSCO_CONFIG_FILE
#

# sga
# python2.7 /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/run_BUSCO.py -i ${SAMPLENAME}_gapfilled.fa --out $SAMPLENAME --lineage_path /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/mammalia_odb9 --mode genome --cpu $CPU --tmp ${PWD}/tmp
python2.7 /mnt/beakerdisk/ANTECHINUS/WGS/BUSCO/busco-master/scripts/run_BUSCO.py -i ${SAMPLENAME}.fa --out $SAMPLENAME --lineage_path /mnt/beakerdisk/ANTECHINUS/WGS/BUSCO/mammalia_odb9 --mode genome --cpu $CPU --tmp ${PWD}/tmp
# Human genome (3.1 Gbp), assessed with 4’104 mammalian BUSCOs: 6 days 15 hours
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# for each
cd $STARTDIR
# 42595A.unique_d10_Q30_genome.fasta
SAMPLENAME=46020A
cd $SAMPLENAME
# need to edit config.ini # ... add BLAST etc. Cmd lines variables below will overrule some of these
# ERROR   No section [busco] found in /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/../config/config.ini. Please make sure both the file and this section exist, see userguide.
#
mkdir -p ./tmp 
export AUGUSTUS_CONFIG_PATH=/bix/augustus-3.2.1/config/ # see https://busco.ezlab.org/v1/files/README.html
export BUSCO_CONFIG_FILE="${PWD}/config.ini"
head $BUSCO_CONFIG_FILE
#

# sga
# python2.7 /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/busco-master/scripts/run_BUSCO.py -i ${SAMPLENAME}_gapfilled.fa --out $SAMPLENAME --lineage_path /mnt/beakerdisk/ANTECHINUS/WGS/short-read/BUSCO/mammalia_odb9 --mode genome --cpu $CPU --tmp ${PWD}/tmp
python2.7 /mnt/beakerdisk/ANTECHINUS/WGS/BUSCO/busco-master/scripts/run_BUSCO.py -i ${SAMPLENAME}.fa --out $SAMPLENAME --lineage_path /mnt/beakerdisk/ANTECHINUS/WGS/BUSCO/mammalia_odb9 --mode genome --cpu $CPU --tmp ${PWD}/tmp
# Human genome (3.1 Gbp), assessed with 4’104 mammalian BUSCOs: 6 days 15 hours
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
