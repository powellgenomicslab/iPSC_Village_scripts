#!/bin/bash
# set -x

# stage="SCTnormalization_no_covatiates"
# stage="SCTnormalization_individual"
# stage="SCTnormalization_CellCycle"
# stage="SCTnormalization_CellCycle_AllTime"
# stage="SCTnormalization_CellCycle_integrated_Time"
# stage="SCTnormalization_CellCycle_integrated_Time_Multi_Resolution"
stage="SCTnormalization_CellCycle_integrated_Time_Clustree"
# stage="SCTnormalization_CellCycle_integrated_Individual"
# stage="SCTnormalization_CellCycle_integrated_Individual_Multi_Resolution"
# stage="SCTnormalization_CellCycle_integrated_Individual_Clustree"


echo $stage

DIR="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
OUT="$DIR/output/Brisbane_Clustering_Testing/$stage/"
ARGS="$OUT/ARGS/"
DATA="$DIR/data/PBMC_scRNA/"
LOG="$OUT/logs/"
PIPELINE="$DIR/scripts/Clustering/Brisbane_Clustering_Testing.sh"


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





if [[ $stage == "SCTnormalization_CellCycle_integrated_Time_Multi_Resolution" ]] || [[ $stage == "SCTnormalization_CellCycle_integrated_Individual_Multi_Resolution" ]]
then
        PE=1
        T=1-16
        MEM=25G
        TEMP=25G
else
        PE=2
        MEM=80G
        TEMP=100G
        T=1-1
fi



qsub -S /bin/bash \
        -q short.q \
        -r yes \
        -l mem_requested=$MEM \
        -l tmp_requested=$TEMP \
        -N Seurat_$stage \
        -t $T \
        -cwd \
        -m e \
        -hold_jid 360112 \
        -M d.neavin@garvan.org.au \
        -j y \
        -e $LOG \
        -o $LOG \
        -V \
        -v DIR=$DIR,ARGfile=$ARGfile,ARGS=$ARGS,stage=$stage,today=$today,OUT=$OUT,DATA=$DATA,LOG=$LOG,PIPELINE=$PIPELINE \
        -C '' $PIPELINE