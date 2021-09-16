#!/bin/bash

PIPELINE=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/Variance_contributions_nb.sh
OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance_contributions_nb/gene_separated/
LOGS=$OUT/logs
mkdir -p $LOGS


qsub -S /bin/bash \
    -q short.q \
    -r yes \
    -t 6979-15292 \
	-tc 250 \
    -l mem_requested=120G \
    -l tmp_requested=120G \
    -N dist_test \
	-hold_jid 181619 \
    -cwd \
    -j y \
    -e $LOGS \
    -o $LOGS \
    -V \
    -v OUT=$OUT \
    -C '' $PIPELINE