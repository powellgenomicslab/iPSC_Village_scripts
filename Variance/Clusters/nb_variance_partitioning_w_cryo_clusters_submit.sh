#!/bin/bash

PIPELINE=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/Clusters/nb_variance_partitioning_w_cryo_clusters.sh
DATA=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/Clusters/nb_variance_partitioning_w_cryo/data/
OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/Clusters/nb_variance_partitioning_w_cryo/gene_separated/
LOGS=$OUT/logs
mkdir -p $LOGS

# number of genes for largest seurat object: 15087


for LOCATION in Brisbane Melbourne Sydney Sydney_Cryopreserved
do
	for CLUSTER in `seq 0 7`
	do
		qsub -S /bin/bash \
			-q short.q \
			-r yes \
			-t 1-15087 \
			-tc 25 \
			-l mem_requested=5G \
			-l tmp_requested=5G \
			-N $LOCATION\_cluster$CLUSTER \
			-cwd \
			-j y \
			-e $LOGS \
			-o $LOGS \
			-V \
			-v OUT=$OUT,LOCATION=$LOCATION,CLUSTER=$CLUSTER,DATA=$DATA \
			-C '' $PIPELINE
	done
done