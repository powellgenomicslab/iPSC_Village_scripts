#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate baseR402

rm /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/hs_err_pid*

Rscript /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/nb_partitioning_village_separate.R $OUT $SGE_TASK_ID $LOCATION $VILLAGE $DATA

rm /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/hs_err_pid*
