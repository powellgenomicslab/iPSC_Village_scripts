#!/bin/bash

PIPELINE=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/nb_partitioning_rep_separate.sh
DATA=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_partitioning_rep_separate/data/
OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/nb_partitioning_rep_separate/gene_separated/
LOGS=$OUT/logs
mkdir -p $LOGS

# number of genes for largest seurat object: 15129


# for REP in Brisbane1 Brisbane2 Brisbane3 Melbourne1 Melbourne2 Melbourne3 Sydney1 Sydney2 Sydney3 Sydney_cryopreserved1 Sydney_cryopreserved2 Sydney_cryopreserved3
for REP in Sydney_cryopreserved1 Sydney_cryopreserved2 Sydney_cryopreserved3 ### Rerun these because coded Time incorrectly
do
	
	echo $REP
	qsub -S /bin/bash \
		-q short.q \
		-r yes \
		-t 1-15129 \
		-tc 75 \
		-hold_jid 190785 \
		-l mem_requested=10G \
		-l tmp_requested=10G \
		-N $REP\_nb \
		-cwd \
		-j y \
		-e $LOGS \
		-o $LOGS \
		-V \
		-v OUT=$OUT,REP=$REP,DATA=$DATA \
		-C '' $PIPELINE
done
