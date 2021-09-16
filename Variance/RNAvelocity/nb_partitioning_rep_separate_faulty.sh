#!/bin/bash

DATA=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning_rep_separate/data/
OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning_rep_separate/gene_separated/
LOGS=$OUT/logs

cd $LOGS

rm $OUT/faulty_runs.txt


for data in `ls $DATA/*.rds`
do
	LOCATION=$(basename $data | cut -d_ -f1-3 | sed 's/_seurat//g')
	echo $LOCATION

	grep 'Model convergence problem' $LOCATION\_nb.o* | cut -d: -f1 >> $OUT/faulty_runs.txt
	grep 'warnings()' $LOCATION\_nb.o* | cut -d: -f1 >> $OUT/faulty_runs.txt
done
