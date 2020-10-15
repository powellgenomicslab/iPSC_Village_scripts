#!/bin/bash
set -x



# stage="CellCycle"
stage="merged"
# stage="mergedRemovedOutliers"

DIR="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
OUT="$DIR/output/ScaterQC/$stage/"
ARGS="$OUT/ARGS/"
DATA="$DIR/data/"
LOG="$OUT/logs/"
EMAIL=d.neavin@garvan.edu.au
PIPELINE="$DIR/scripts/Scater/Scater.sh"

today=`date +%Y-%m-%d.%H:%M:%S`
ARGfile="$ARGS/arguments$stage$today.txt"

mkdir -p $ARGS
mkdir -p $LOG
mkdir -p $OUT

dir -d $DATA/200128_A00152_0196_BH3HNFDSXY/GE/DRENEA*/outs/filtered_feature_bc_matrix/ > $DATA/SampleDIRS.txt

echo $DATA > $ARGfile
echo $OUT >> $ARGfile
echo $stage >> $ARGfile

export PATH=$PATH

cd $OUT

if [[ $stage == "CellCycle" ]]
then 
    MEM=30G
    NUMSAMPS=$(cat $DATA/SampleDIRS.txt | wc -l)
    TC=$NUMSAMPS
elif [[ $stage == "merged" ]] || [[ $stage == "mergedRemovedOutliers" ]]
then
    NUMSAMPS=1
    TC=1
    MEM=200G
fi
### Note: use 2-$NUMSAMPS below because Pool 10 is unwanted and is number 1
qsub -S /bin/bash \
        -q short.q \
        -r yes \
        -t 1-$NUMSAMPS \
        -tc $TC \
        -l mem_requested=$MEM \
        -l tmp_requested=$MEM \
        -N Scater_$stage \
        -cwd \
        -m e \
        -M d.neavin@garvan.org.au \
        -j y \
        -e $LOG \
        -o $LOG \
        -V \
        -v DIR=$DIR,ARGfile=$ARGfile,ARGS=$ARGS,stage=$stage,today=$today \
        -C '' $PIPELINE

#        
#