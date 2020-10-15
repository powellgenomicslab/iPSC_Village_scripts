#!/bin/bash

DIR="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village"
DATA="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/200128_A00152_0196_BH3HNFDSXY/GE/"
OUT="$DIR/output/scSplit"
LOG="$OUT/logs"
PIPELINE="/$DIR/scripts/scSplit/scSplit.sh"

# stage="prepare_files"
# stage="freebayes"
# stage="AlleleMatrices"
stage="Demultiplex"
# stage="scSplitGenotypes"
echo $stage

if [[ $stage == "prepare_files" ]]
then
    T=1-6
    PE=40
    TC=6
    MEM=1G
elif [[ $stage == "freebayes" ]]
then
    PE=40
    T=1-10
    TC=10
    MEM=1G
elif [[ $stage == "AlleleMatrices" ]]
then
    PE=5
    T=1-6
    TC=10
    MEM=50G
elif [[ $stage == "Demultiplex" ]]
then
    PE=10
    T=1-10
    TC=10
    MEM=50G
elif [[ $stage == "scSplitGenotypes" ]]
then
    PE=1
    T=8
    TC=1-10
    MEM=400G
fi

mkdir -p $OUT
mkdir -p $LOG

qsub -S /bin/bash \
    -q short.q \
    -r yes \
    -t $T \
    -tc $TC \
    -l mem_requested=$MEM \
    -l tmp_requested=$MEM \
    -pe smp $PE \
    -N scSplit_$stage \
    -cwd \
    -m e \
    -M d.neavin@garvan.org.au \
    -j y \
    -e $LOG \
    -o $LOG \
    -V \
    -hold_jid 337826 \
    -v DATA=$DATA,DIR=$DIR,OUT=$OUT,FASTA=$FASTA,today=$today,stage=$stage,PATH=$PATH,stage=$stage,T=$PE \
    -C '' $PIPELINE