#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate baseR402


Rscript /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/RNAvelocity/nb_partitioning_rep_separate.R $OUT $SGE_TASK_ID $REP $DATA

# rm /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/hs_err_pid*
