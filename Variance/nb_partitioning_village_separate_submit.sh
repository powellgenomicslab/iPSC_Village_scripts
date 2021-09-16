#!/bin/bash

PIPELINE=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/nb_partitioning_village_separate.sh
DATA=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_partitioning_village_separate/data/
OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_partitioning_village_separate/gene_separated/
LOGS=$OUT/logs
mkdir -p $LOGS

# number of genes for largest seurat object: 16258


for LOCATION in Brisbane Melbourne Sydney_Fresh Sydney_Cryopreserved
do
	echo $LOCATION
	for VILLAGE in Baseline Village
	do
		qsub -S /bin/bash \
			-q short.q \
			-r yes \
			-t 1-16258 \
			-tc 25 \
			-l mem_requested=10G \
			-l tmp_requested=10G \
			-N $LOCATION\_$VILLAGE\_nb \
			-cwd \
			-j y \
			-e $LOGS \
			-o $LOGS \
			-V \
			-v OUT=$OUT,LOCATION=$LOCATION,VILLAGE=$VILLAGE,DATA=$DATA \
			-C '' $PIPELINE
	done
done