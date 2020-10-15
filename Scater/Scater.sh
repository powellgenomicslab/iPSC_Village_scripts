#!/bin/bash
RSCRIPT="$DIR/scripts/Scater/Scater.R"

cat $ARGfile > $ARGS/arguments$stage$today.$SGE_TASK_ID.txt
echo $SGE_TASK_ID >> $ARGS/arguments$stage$today.$SGE_TASK_ID.txt

eval "$(conda shell.bash hook)"
conda activate generalR
Rscript $RSCRIPT $ARGS/arguments$stage$today.$SGE_TASK_ID.txt
conda deactivate