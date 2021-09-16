#!/bin/bash


SCRIPT="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/RNAvelocity/velocyto.sh"
TENxDIR="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE"
OUT="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/preprocess"
GTF="/directflow/SCCGGroupShare/projects/DrewNeavin/References/GenCode/GRCh38/gencode.v38.annotation.gtf"
LOG=$OUT/logs

T=32


# for pool in DRENEA_{1..6} Village_B_1_week Village_A_Baseline
for pool in DRENEA_1
do
	qsub -S /bin/bash \
        -q short.q \
        -r yes \
        -l mem_requested=4G \
        -l tmp_requested=4G \
		-pe smp $T \
        -N $pool\_velocyto \
        -cwd \
        -m e \
        -M d.neavin@garvan.org.au \
        -j y \
        -e $LOG \
        -o $LOG \
        -V \
        -v TENxDIR=$TENxDIR,OUT=$OUT,pool=$pool,GTF=$GTF,T=$T \
        -C '' $SCRIPT
done