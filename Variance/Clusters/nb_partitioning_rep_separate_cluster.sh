#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate baseR402


Rscript /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/Clusters/nb_partitioning_rep_separate_cluster.R $OUT $SGE_TASK_ID $REP $CLUSTER $DATA

rm /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/Clusters/hs_err_pid*
