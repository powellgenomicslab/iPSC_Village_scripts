#!/bin/bash

eval "$(conda shell.bash hook)"
conda activate baseR402

Rscript /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/Distribution_tests_parallel.R $OUT $SGE_TASK_ID