#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate baseR402

Rscript /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/nb_variance_partitioning.R $OUT $SGE_TASK_ID