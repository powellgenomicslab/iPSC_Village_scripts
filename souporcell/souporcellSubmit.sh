#!/bin/bash
set -x

stage="run_souporcell"
# stage="fix_output_files"

DIR="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
DATA="$DIR/data/200128_A00152_0196_BH3HNFDSXY/GE/"
OUT="$DIR/output/souporcell"
LOG="$OUT/logs"
PIPELINE=$DIR/scripts/souporcell/souporcell.sh

mkdir -p $OUT
mkdir -p $LOG



if [[ $stage == "run_souporcell" ]]
then
    T=1-6
    PE=20
    MEM=10G
    Q=short.q
    TC=3
elif [[ $stage == "fix_output_files" ]]
then
    T=1-6
    PE=1
    MEM=10G
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
    -N souporcell_$stage \
    -cwd \
    -m e \
    -M d.neavin@garvan.org.au \
    -hold_jid 308370 \
    -j y \
    -e $LOG \
    -o $LOG \
    -v T=$PE,DATA=$DATA,DIR=$DIR,OUT=$OUT,today=$today,stage=$stage,ARGS=$ARGS \
    -C '' $PIPELINE

