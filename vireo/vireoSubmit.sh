#!/bin/bash


# FASTA=/share/ScratchGeneral/drenea/References/GenCode/GRCh38/gencode.v30.transcripts.fa
DIR="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
DATA="$DIR/data/200128_A00152_0196_BH3HNFDSXY/GE/"
OUT="$DIR/output/vireo"
LOG="$OUT/logs"
PIPELINE=$DIR/scripts/vireo/vireo.sh
ARGS="$OUT/ARGS/" 

# stage="prepBAMS" 
# stage="pileup" 
# stage="vireoPrepareImputed"
stage="vireoGenotyped"
# stage="vireoUngenotyped"
 
mkdir -p $ARGS
mkdir -p $OUT
mkdir -p $LOG

today=`date +%Y-%m-%d.%H:%M:%S`
ARGfile="$ARGS/arguments$today.txt"


if [[ $stage == "prepBAMS" ]] || [[ $stage == "pileup" ]] || [[ $stage == "vireoPrepareImputed" ]]
then
    NUMSAMPS=10
    T=1-6
    PE=7
    MEM=10G
    Q=short.q
    TC=6
elif [[ $stage == "vireoGenotyped" ]]  || [[ $stage == "vireoUngenotyped" ]]
then
    T=1-6
    NUMSAMPS=10
    PE=4
    MEM=100G
    Q=short.q
    TC=10
fi

qsub -S /bin/bash \
    -q $Q \
    -r yes \
    -t $T \
    -tc $TC \
    -l mem_requested=$MEM \
    -l tmp_requested=$MEM \
    -pe smp $PE \
    -N vireo_$stage \
    -cwd \
    -m e \
    -M d.neavin@garvan.org.au \
    -hold_jid 308370 \
    -j y \
    -e $LOG \
    -o $LOG \
    -v PE=$PE,DATA=$DATA,DIR=$DIR,OUT=$OUT,today=$today,stage=$stage,ARGS=$ARGS \
    -C '' $PIPELINE

