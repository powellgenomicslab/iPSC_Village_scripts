#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate baseR402

Rscript /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/Variance_contributions_nb.R $OUT $SGE_TASK_ID