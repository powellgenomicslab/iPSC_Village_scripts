#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate baseR402

rm /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/Clusters/hs_err_pid*

Rscript /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/Clusters/nb_variance_partitioning_w_cryo_clusters.R $OUT $SGE_TASK_ID $LOCATION $CLUSTER $DATA

rm /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/Clusters/hs_err_pid*
