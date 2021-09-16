#!/bin/bash


SCRIPT="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/RNAvelocity/bam2fastq.sh"
TENxDIR="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE"
OUT="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/preprocess"
LOG=$OUT/logs

mkdir -p $LOG


# for pool in DRENEA_{2..6} Village_B_1_week Village_A_Baseline
for pool in DRENEA_1
do
	qsub -S /bin/bash \
        -q short.q \
        -r yes \
        -l mem_requested=50G \
        -l tmp_requested=50G \
        -N $pool\_bam2fastq \
        -cwd \
        -m e \
        -M d.neavin@garvan.org.au \
        -j y \
        -e $LOG \
        -o $LOG \
        -V \
        -v TENxDIR=$TENxDIR,OUT=$OUT,pool=$pool \
        -C '' $SCRIPT
done