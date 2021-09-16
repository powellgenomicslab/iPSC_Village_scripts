#!/bin/bash

DATA=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_partitioning_cryo/data/
OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_variance_partitioning_w_cryo/gene_separated/
LOGS=$OUT/logs

cd $LOGS

rm $OUT/faulty_runs.txt

for data in `ls $DATA/*.rds`
do
	LOCATION=$(basename $data | cut -d_ -f1-2 | sed 's/_SCT//g')
	echo $LOCATION

	grep 'Model convergence problem' $LOCATION\_nb.o* | cut -d: -f1 >> $OUT/faulty_runs.txt
	grep 'warnings()' $LOCATION\_nb.o* | cut -d: -f1 >> $OUT/faulty_runs.txt
done
