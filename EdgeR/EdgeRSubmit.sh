#!/bin/bash

# stage="PrePost"
# stage="PrePost_DE"
# stage="PrePost_All"
# stage="PrePost_All_DE"
# stage="PrePost_CellLines"
stage="PrePost_CellLines_DE"

DIR="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
DATA="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Seurat/Harmony_DayPatient/"
OUT="$DIR/output/EdgeR/$stage/"
ARGS="$OUT/ARGS/"
LOG="$OUT/logs/"
PIPELINE="$DIR/scripts/EdgeR/EdgeR.sh"

today=`date +%Y-%m-%d.%H:%M:%S`
ARGfile="$ARGS/arguments$stage$today.txt"

mkdir -p $ARGS
mkdir -p $LOG
mkdir -p $OUT

echo $stage > $ARGfile
echo $OUT >> $ARGfile
echo $DIR >> $ARGfile

export PATH=$PATH

cd $OUT

 ### For new stages ###
if [[ $stage == "PrePost" ]] || [[ $stage == "PrePost_CellLines" ]]
then
    MEM=30G
    PE=2
    T=1-3
elif [[ $stage == "PrePost_DE" ]] || [[ $stage == "PrePost_CellLines_DE" ]]
then
    MEM=1G
    PE=100
    T=1-3
elif [[ $stage == "PrePost_All" ]]
then
    MEM=30G
    PE=2
    T=1-1
elif [[ $stage == "PrePost_All_DE" ]]
then
    MEM=1G
    PE=100
    T=1-1
fi



qsub -S /bin/bash \
        -q short.q \
        -r yes \
        -l mem_requested=$MEM \
        -l tmp_requested=$MEM \
        -l h=!*beefy-4-4* \
        -l h=!*epsilon* \
        -N EdgeR_$stage \
        -pe smp $PE \
        -t $T \
        -cwd \
        -m e \
        -hold_jid 312760 \
        -M d.neavin@garvan.org.au \
        -j y \
        -e $LOG \
        -o $LOG \
        -V \
        -v DIR=$DIR,ARGfile=$ARGfile,ARGS=$ARGS,stage=$stage,today=$today,OUT=$OUT,DATA=$DATA,LOG=$LOG,PIPELINE=$PIPELINE \
        -C '' $PIPELINE