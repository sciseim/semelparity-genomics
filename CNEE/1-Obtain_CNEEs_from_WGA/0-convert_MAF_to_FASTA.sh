#!/usr/bin/env bash

# convert MAF to FASTA
for MAF in *.maf; do cat $MAF | awk '/^s/{print ">" $2 "\n" $7}' >>$MAF.fa;done
