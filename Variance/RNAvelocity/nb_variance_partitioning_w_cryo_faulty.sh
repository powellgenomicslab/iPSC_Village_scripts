#!/bin/bash

DATA=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning/data/
OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning/gene_separated/
LOGS=$OUT/logs

cd $LOGS

rm $OUT/faulty_runs.txt

for data in `ls $DATA/*.rds`
do
	LOCATION=$(basename $data | cut -d_ -f1-3 | sed 's/_SCT//g')
	echo $LOCATION

	grep 'Model convergence problem' $LOCATION\_nb.o* | cut -d: -f1 >> $OUT/faulty_runs.txt
	grep 'warnings()' $LOCATION\_nb.o* | cut -d: -f1 >> $OUT/faulty_runs.txt
done
