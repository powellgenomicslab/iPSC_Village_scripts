#!/bin/bash

PIPELINE=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/nb_variance_partitioning.sh
OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_variance_partitioning/gene_separated/
LOGS=$OUT/logs
mkdir -p $LOGS


qsub -S /bin/bash \
    -q short.q \
    -r yes \
    -t 1-15292 \
	-tc 150 \
    -l mem_requested=120G \
    -l tmp_requested=120G \
    -N nb_dist_test \
    -cwd \
    -j y \
    -e $LOGS \
    -o $LOGS \
    -V \
    -v OUT=$OUT \
    -C '' $PIPELINE