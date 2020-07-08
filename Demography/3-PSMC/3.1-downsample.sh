# #!/usr/bin/env bash

# downsample A. flavipes samples to 30X (to match the other samples)
# samtools view -h -s 0.3 in.bam >10x.bam # to obtain 10X if starting was 33X
# AdamAnt 44x 30/44 = 0.6818182
# Ant2 41x 30/41=0.7317073

# thus ... 
cd /mnt/beakerdisk/ANTECHINUS/PSMC/psmc/TasDevNoChrX/downsampled_BAMs
OLDSAMTOOLSPATH=/mnt/beakerdisk/ANTECHINUS/WGS/short-read/reference-based/TassieDevil/RefGenBin/samtools-0.1.19
time ${OLDSAMTOOLSPATH}/samtools view -s 0.68 /mnt/beakerdisk/ANTECHINUS/PSMC/pseudogenomes-no-chrX/AdamAnt/AdamAnt.aln.sorted.unique.bam >AdamAnt_30x.bam && \
time ${OLDSAMTOOLSPATH}/samtools view -s 0.68 /mnt/beakerdisk/ANTECHINUS/PSMC/pseudogenomes-no-chrX/Ant2/Ant2.aln.sorted.unique.bam >Ant2_30x.bam 