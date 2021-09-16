#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate baseR402

rm /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/hs_err_pid*

Rscript /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/nb_variance_partitioning_w_cryo.R $OUT $SGE_TASK_ID $LOCATION $DATA

rm /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/hs_err_pid*
