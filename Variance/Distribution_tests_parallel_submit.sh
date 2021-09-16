#!/bin/bash

##### Do this before running for first time

PIPELINE=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/Distribution_tests_parallel.sh
OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Distribution_tests_parallel/gene_separated/
LOGS=$OUT/logs
mkdir -p $LOGS


qsub -S /bin/bash \
    -q short.q \
    -r yes \
    -t 1-15292 \
	-tc 250 \
    -l mem_requested=40G \
    -l tmp_requested=40G \
    -pe smp 1 \
    -N dist_test \
    -cwd \
    -j y \
    -e $LOGS \
    -o $LOGS \
    -V \
    -v OUT=$OUT \
    -C '' $PIPELINE