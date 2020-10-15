#!/bin/bash
RSCRIPT="$DIR/scripts/Seurat/Seurat_w_Hash.R"

echo "The stage is: $stage"
echo "The main directory is: $DIR"
echo "The out directory is: $OUT"
echo "The argument directory is: $ARGS"
echo "The data directory is: $DATA"
echo "The log directory is: $LOG"
echo "The R script being used is: $PIPELINE"
echo "The date is: $today"
echo "The arguments file is: $ARGfile"

cat $ARGfile > $ARGS/arguments$stage$today.$SGE_TASK_ID.txt
echo $SGE_TASK_ID >> $ARGS/arguments$stage$today.$SGE_TASK_ID.txt

##### If running on wolfpack 7 #####
if [[ $stage == "CellCycle" ]]
then
    eval "$(conda shell.bash hook)"
    conda activate harmony
else
    eval "$(conda shell.bash hook)"
    conda activate harmony
fi

Rscript $RSCRIPT $ARGS/arguments$stage$today.$SGE_TASK_ID.txt