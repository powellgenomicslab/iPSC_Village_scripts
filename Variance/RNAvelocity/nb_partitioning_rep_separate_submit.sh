#!/bin/bash

PIPELINE=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/RNAvelocity/nb_partitioning_rep_separate.sh
DATA=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning_rep_separate/data/
OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/nb_partitioning_rep_separate/gene_separated/
LOGS=$OUT/logs
mkdir -p $LOGS

# number of genes for largest seurat object: 15106


for FILE in `ls $DATA/*.rds`
do
	REP=$(basename $FILE | cut -d_ -f1-3 | sed 's/_seurat//g')
	
	echo $REP

	qsub -S /bin/bash \
		-q short.q \
		-r yes \
		-t 1-15106 \
		-tc 50 \
		-hold_jid 292644 \
		-l mem_requested=2G \
		-l tmp_requested=2G \
		-N $REP\_nb \
		-cwd \
		-j y \
		-e $LOGS \
		-o $LOGS \
		-V \
		-v OUT=$OUT,REP=$REP,DATA=$DATA \
		-C '' $PIPELINE
done


