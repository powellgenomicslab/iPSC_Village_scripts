#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate baseR402

OUTDIR=$OUT/$LOCATION/

mkdir $OUTDIR

rm /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/RNAvelocity/hs_err_pid*

Rscript /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/RNAvelocity/nb_variance_partitioning_w_cryo_latent.R $OUTDIR $SGE_TASK_ID $LOCATION $DATA

rm /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/RNAvelocity/hs_err_pid*
