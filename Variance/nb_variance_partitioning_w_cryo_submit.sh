#!/bin/bash

PIPELINE=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/nb_variance_partitioning_w_cryo.sh
DATA=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_partitioning_cryo/data/
OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_variance_partitioning_w_cryo/gene_separated/
LOGS=$OUT/logs
mkdir -p $LOGS

# number of genes for largest seurat object: 15127


# for LOCATION in Brisbane Melbourne Sydney Sydney_Cryopreserved
for LOCATION in Sydney_Cryopreserved ### Redo because accidently didn't have the Time correctly coded
do
	qsub -S /bin/bash \
		-q short.q \
		-r yes \
		-t 1-15127 \
		-tc 200 \
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