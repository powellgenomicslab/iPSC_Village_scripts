#!/bin/bash

PIPELINE=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/RNAvelocity/nb_variance_partitioning_w_cryo.sh
DATA=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning/data/
OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning/gene_separated/
LOGS=$OUT/logs
mkdir -p $LOGS

# number of genes for largest seurat object: 15119


for FILE in `ls $DATA/*.rds`
do
	LOCATION=$(basename $FILE | cut -d_ -f1-3 | sed 's/_SCT//g')
	
	echo $LOCATION

	
	qsub -S /bin/bash \
		-q short.q \
		-r yes \
		-t 1-15119 \
		-tc 200 \
		-l mem_requested=4G \
		-l tmp_requested=4G \
		-N $LOCATION\_nb \
		-cwd \
		-j y \
		-e $LOGS \
		-o $LOGS \
		-V \
		-v OUT=$OUT,LOCATION=$LOCATION,DATA=$DATA \
		-C '' $PIPELINE
done


