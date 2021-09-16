#!/bin/bash


SCRIPT="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/RNAvelocity/salmon_alevin.sh"
OUT="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/preprocess"
LOG=$OUT/logs

mkdir -p $LOG

T=10

for pool in DRENEA_{1..6} Village_B_1_week Village_A_Baseline
do
	qsub -S /bin/bash \
        -q short.q \
        -r yes \
        -l mem_requested=4G \
        -l tmp_requested=4G \
        -N $pool\_alevin \
		-pe smp $T \
        -cwd \
        -m e \
        -M d.neavin@garvan.org.au \
        -j y \
        -e $LOG \
        -o $LOG \
        -V \
        -v OUT=$OUT,pool=$pool,T=$T \
        -C '' $SCRIPT
done