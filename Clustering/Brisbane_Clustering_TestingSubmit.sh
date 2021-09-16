#!/bin/bash
# set -x

# stage="SCTnormalization_no_covatiates"
# stage="SCTnormalization_individual"
# stage="SCTnormalization_CellCycle"
# stage="SCTnormalization_CellCycle_AllTime"
# stage="SCTnormalization_CellCycle_integrated_Time"
# stage="SCTnormalization_CellCycle_integrated_Time_Multi_Resolution"
# stage="SCTnormalization_CellCycle_integrated_Time_Clustree"
# stage="SCTnormalization_CellCycle_integrated_Individual"
# stage="SCTnormalization_CellCycle_integrated_Individual_Multi_Resolution"
# stage="SCTnormalization_CellCycle_integrated_Individual_Clustree"
# stage="SCTnormalization_Baseline_CellCycle_integrated_Individual"
# stage="SCTnormalization_Baseline_CellCycle_integrated_Individual_Multi_Resolution"
# stage="SCTnormalization_Baseline_CellCycle_integrated_Individual_Clustree"
# stage="SCTnormalization_Baseline_CellCycle_integrated_Individual_DE"
# stage="Make_Final_Brisbane_Seurat"
# stage="Make_Final_Brisbane_Baseline_Seurat"
# stage="SCTnormalization_CellCycle_Pool_integrated_Time"
# stage="SCTnormalization_Baseline_CellCycle_Pool_integrated_Time"
# stage="SCTnormalization_Baseline_CellCycle_Pool_mt_integrated_Time"
# stage="SCTnormalization_Baseline_CellCycle_Pool_mt_integrated_Time_Multi_Resolution"
# stage="SCTnormalization_Baseline_CellCycle_Pool_mt_integrated_Time_Clustree"
# stage="SCTnormalization_Baseline_CellCycle_Pool_mt_integrated_Time_DE"
stage="Make_Brisbane_Baseline_Seurat_covs_res_0.13"
# stage="Make_Brisbane_Baseline_Seurat_covs_res_0.2"

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





if [[ $stage == "SCTnormalization_CellCycle_integrated_Time_Multi_Resolution" ]] || [[ $stage == "SCTnormalization_CellCycle_integrated_Individual_Multi_Resolution" ]] || [[ $stage == "SCTnormalization_Baseline_CellCycle_integrated_Individual_Multi_Resolution" ]] || [[ $stage == "SCTnormalization_Baseline_CellCycle_Pool_mt_integrated_Time_Multi_Resolution" ]]
then
        PE=1
        T=1-16
        T=27-36
        MEM=25G
        TEMP=25G
else
        PE=2
        MEM=80G
        TEMP=100G
        T=1-1
fi

MEM=5G
TEMP=5G


qsub -S /bin/bash \
        -r yes \
        -l mem_requested=$MEM \
        -l tmp_requested=$TEMP \
        -N Seurat_$stage \
        -t $T \
        -cwd \
        -m e \
        -M d.neavin@garvan.org.au \
        -j y \
        -e $LOG \
        -o $LOG \
		-q long.q \
        -V \
        -v DIR=$DIR,ARGfile=$ARGfile,ARGS=$ARGS,stage=$stage,today=$today,OUT=$OUT,DATA=$DATA,LOG=$LOG,PIPELINE=$PIPELINE \
        -C '' $PIPELINE