#!/bin/bash

. /etc/profile.d/modules.sh # This is a line that lets you load modules in this script and will include them in the qsub submission

# stage="pileup"
# stage="freemuxlet"
stage="demuxlet"

DIR="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
DATA="$DIR/data/200128_A00152_0196_BH3HNFDSXY/GE/"
OUT="$DIR/output/popscle/$stage/"
LOGS="$OUT/logs/"
PIPELINE=$DIR/scripts/popscle/freemuxlet.sh

mkdir -p $LOGS
mkdir -p $OUT



export $PATH

if [[ $stage == "pileup" ]]
then
    MEM=25G
    PE=20
    TC=10
    T=1-6
elif [[ $stage == "demuxlet" ]]
then
    MEM=10G
    TC=10
    T=1-6
    PE=2
elif [[ $stage == "freemuxlet" ]]
then
    MEM=25G
    TC=10
    T=1-6
    PE=2
fi

qsub -S /bin/bash \
    -q short.q \
    -r yes \
    -l mem_requested=$MEM \
    -l tmp_requested=$MEM \
    -hold_jid 334378 \
    -N Freemuxlet_$stage \
    -M d.neavin@garvan.org.au \
    -t $T \
    -tc $TC \
    -pe smp $PE \
    -j y \
    -e $LOGS \
    -o $LOGS \
    -m e \
    -M d.neavin@garvan.org.au \
    -v DIR=$DIR,OUT=$OUT,stage=$stage,PATH=$PATH,DATA=$DATA \
    -V \
    -C '' $PIPELINE
