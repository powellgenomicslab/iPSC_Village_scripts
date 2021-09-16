#!/bin/bash

PIPELINE=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/Clusters/nb_partitioning_rep_separate_cluster.sh
DATA=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/Clusters/nb_partitioning_rep_separate/data/
OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/Clusters/nb_partitioning_rep_separate/gene_separated/
LOGS=$OUT/logs
mkdir -p $LOGS

# number of genes for largest seurat object: 15077


for REP in Brisbane1 Brisbane2 Brisbane3 Melbourne1 Melbourne2 Melbourne3 Sydney1 Sydney2 Sydney3 Sydney_cryopreserved1 Sydney_cryopreserved2 Sydney_cryopreserved3
do
	for CLUSTER in `seq 0 7`
	do
	echo $REP
	qsub -S /bin/bash \
		-q short.q \
		-r yes \
		-t 1-15077 \
		-tc 20 \
		-l mem_requested=2G \
		-l tmp_requested=2G \
		-N $REP\_cluster$CLUSTER \
		-cwd \
		-j y \
		-e $LOGS \
		-o $LOGS \
		-V \
		-v OUT=$OUT,REP=$REP,CLUSTER=$CLUSTER,DATA=$DATA \
		-C '' $PIPELINE
	done
done
