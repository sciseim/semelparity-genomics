#!/usr/bin/env bash
PSMCPATH=/mnt/beakerdisk/ANTECHINUS/WGS/short-read/reference-based/TassieDevil/RefGenBin/psmc-master

SAMPLENAME=BD-17-5A
cd /mnt/beakerdisk/ANTECHINUS/PSMC/psmc/TasDevNoChrX/${SAMPLENAME}

# 100 bootstraps
# bootstrapping in parallel
nohup parallel -j15 --nice 19 "${PSMCPATH}/psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" ${SAMPLENAME}_split.psmcfa -o round-{}.psmc" :::: <(seq 100) &
# 10 jobs at the same, so -- with 20 cores -- I can run two parallel sessions (e.g. two species)
