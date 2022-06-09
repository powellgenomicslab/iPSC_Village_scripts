#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate baseR402


Rscript /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/variance_partition_post_review.R $OUT_ICC $OUT_MODEL $OUT_RESID $SGE_TASK_ID