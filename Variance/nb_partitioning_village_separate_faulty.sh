#!/bin/bash

DATA=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_partitioning_village_separate/data
OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_partitioning_village_separate/gene_separated/
LOGS=$OUT/logs

cd $LOGS

rm $OUT/faulty_runs.txt

for data in `ls $DATA/*.rds`
do
	LOCATION=$(basename $data | sed 's/_seurat.rds//g' | sed 's/Village_//g' | sed 's/Baseline_//g')
	TIME=$(basename $data | sed 's/_seurat.rds//g' | sed "s/_$LOCATION//g")
	echo $LOCATION
	echo $TIME

	grep 'Model convergence problem' $LOCATION\_$TIME\_nb.o* | cut -d: -f1 >> $OUT/faulty_runs.txt
	grep 'warnings()' $LOCATION\_$TIME\_nb.o* | cut -d: -f1 >> $OUT/faulty_runs.txt
done
