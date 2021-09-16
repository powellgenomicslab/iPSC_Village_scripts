#!/bin/bash

PIPELINE=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/RNAvelocity/nb_variance_partitioning_w_cryo_latent.sh
DATA=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning_w_latent/data/
OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning_w_latent/gene_separated/
LOGS=$OUT/logs
mkdir -p $LOGS

# number of genes for largest seurat object: 15129


for FILE in `ls $DATA/*.rds`
do
	LOCATION=$(basename $FILE | cut -d_ -f1-2 | sed 's/_SCT//g')
	
	echo $LOCATION

	
	qsub -S /bin/bash \
		-q short.q \
		-r yes \
		-t 1-15129 \
		-tc 400 \
		-l mem_requested=15G \
		-l tmp_requested=15G \
		-N $LOCATION\_nb \
		-cwd \
		-j y \
		-e $LOGS \
		-o $LOGS \
		-V \
		-v OUT=$OUT,LOCATION=$LOCATION,DATA=$DATA \
		-C '' $PIPELINE
done