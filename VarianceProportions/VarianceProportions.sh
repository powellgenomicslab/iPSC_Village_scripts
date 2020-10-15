RSCRIPT="$DIR/scripts/VarianceProportions/VarianceProportions.R"

cat $ARGfile > $ARGS/arguments$stage$today.$SGE_TASK_ID.txt
echo $SGE_TASK_ID >> $ARGS/arguments$stage$today.$SGE_TASK_ID.txt

echo "The stage is: $stage"
echo "The main directory is: $DIR"
echo "The out directory is: $OUT"
echo "The argument directory is: $ARGS"
echo "The data directory is: $DATA"
echo "The log directory is: $LOG"
echo "The R script being used is: $PIPELINE"
echo "The date is: $today"
echo "The arguments file is: $ARGfile"


##### If running on wolfpack 7 #####
eval "$(conda shell.bash hook)"
conda activate generalR

/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/harmony/lib/R/bin/Rscript $RSCRIPT $ARGS/arguments$stage$today.$SGE_TASK_ID.txt