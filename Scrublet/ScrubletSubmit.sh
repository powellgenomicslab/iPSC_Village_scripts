#!/bin/bash
set -x

echo "Loading modules"
. /etc/profile.d/modules.sh
module load /share/ClusterShare/Modules/modulefiles/contrib/briglo/miniconda/3

echo "Starting Script"
DIR="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village"

name="default0.85"
# name="Var0.8"
# name="Var0.90"
# name="Var0.95"

echo "Assigning Variables"
OUT="$DIR/output/Scrublet/$name"
PIPELINE="$DIR/scripts/Scrublet/Scrublet.sh"
SAMPLEINFO="$DIR/data/SampleInfo.txt"

NUMSAMPS=$(cat $SAMPLEINFO | wc -l)

echo "Making Directories"
mkdir -p $OUT
mkdir -p $OUT/logs

cd $OUT

if [[ $name == "default0.85" ]]
then
    VAR=0.85
elif [[ $name == "Var0.8" ]]
then
    VAR=0.8
elif [[ $name == "Var0.90" ]]
then
    VAR=0.9
elif [[ $name == "Var0.95" ]]
then
    VAR=0.95
fi

NUMSAMPS=1

echo "Starting qsub"
qsub -S /bin/bash \
    -cwd \
    -t 1-$NUMSAMPS \
    -N array \
    -q short.q \
    -l mem_requested=20G \
    -l tmp_requested=20G \
    -e $OUT/logs \
    -o $OUT/logs \
    -r yes \
    -j y \
    -M d.neavin@garvan.org.au \
    -v OUT=$OUT,SAMPLEINFO=$SAMPLEINFO,DIR=$DIR,VAR=$VAR \
    -V \
    -C '' $PIPELINE