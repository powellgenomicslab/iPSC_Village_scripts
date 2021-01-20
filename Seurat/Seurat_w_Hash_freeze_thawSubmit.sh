#!/bin/bash
set -x

# stage="Hashtags"
# stage="Demultiplex_Doublet"
# stage="QC"
stage="CellCycle"
# stage="Old_scRNAseq_data"
# stage="SCTnormalization_no_covatiates"
# stage="SCTnormalization_PoolRegressed"
# stage="SCTnormalization_PoolDayRegressed"
# stage="SCTnormalization_SiteRegressed"
# stage="SCTnormalization_SiteMtRbRegressed"
# stage="Harmony_Pool"
# stage="Harmony_Site"
# stage="Harmony_Site_Regress_MtRb"
# stage="Harmony_Day"
# stage="Harmony_DayPatient"
# stage="Integration_DayPatient"
# stage="Integration_Day"
# stage="Integration_Site"
# stage="Integration_Site_Regress_MtRb"
# stage="Integration_Day_then_Individual"
# stage="Multi_Resolution_no_covatiates"
# stage="Multi_Resolution_PoolRegressed"
# stage="Multi_Resolution_PoolDayRegressed"
# stage="Multi_Resolution_Harmony_Pool"
# stage="Multi_Resolution_Harmony_Day"
# stage="Multi_Resolution_Harmony_DayPatient"
# stage="Multi_Resolution_Integration_DayPatient"
# stage="Multi_Resolution_Integration_Day_then_Individual"
# stage="Clustree_no_covatiates"
# stage="Clustree_PoolRegressed"
# stage="Clustree_PoolDayRegressed"
# stage="Clustree_Harmony_Pool"
# stage="Clustree_Harmony_Day"
# stage="Clustree_Harmony_DayPatient"
# stage="Clustree_Integration_DayPatient"
# stage="Clustree_Integration_Day_then_Individual"
# stage="QCfigurew_no_covatiates"
# stage="QCfigurew_PoolRegressed"
# stage="QCfigurew_PoolDayRegressed"
# stage="QCfigures_Harmony_Pool"
# stage="QCfigures_Harmony_Day"
# stage="QCfigures_Harmony_DayPatient"
# stage="QCfigures_Integration_DayPatient"
# stage="QCfigures_Integration_Day_then_Individual"

DIR="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
OUT="$DIR/output/Seurat_w_Hash_freeze_thaw/$stage/"
ARGS="$OUT/ARGS/"
DATA="$DIR/data/PBMC_scRNA/"
LOG="$OUT/logs/"
PIPELINE="$DIR/scripts/Seurat/Seurat_w_Hash_freeze_thaw.sh"

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

if [[ $stage == "Hashtags" ]] || [[ $stage == "QC" ]] || [[ $stage == "CellCycle" ]] || [[ $stage == "Old_scRNAseq_data" ]] || [[ $stage == "SCTnormalization_no_covatiates" ]] || [[ $stage == "SCTnormalization_SiteRegressed" ]] || [[ $stage == "SCTnormalization_PoolRegressed" ]] || [[ $stage == "Clustree_no_covatiates" ]] || [[ $stage == "Clustree_PoolRegressed" ]] || [[ $stage == "QCfigurew_no_covatiates" ]] || [[ $stage == "QCfigurew_PoolRegressed" ]] || [[ $stage == "SCTnormalization_PoolDayRegressed" ]] || [[ $stage == "Clustree_PoolDayRegressed" ]] || [[ $stage == "QCfigurew_PoolDayRegressed" ]] || [[ $stage == "Harmony_Pool" ]] || [[ $stage == "Harmony_Site" ]] || [[ $stage == "Harmony_Site_Regress_MtRb" ]] || [[ $stage == "Harmony_Day" ]] || [[ $stage == "Clustree_Harmony_Pool" ]] || [[ $stage == "Clustree_Harmony_Day" ]] || [[ $stage == "QCfigures_Harmony_Pool" ]] || [[ $stage == "QCfigures_Harmony_Day" ]] || [[ $stage == "Harmony_DayPatient" ]] || [[ $stage == "Clustree_Harmony_DayPatient" ]] || [[ $stage == "QCfigures_Harmony_DayPatient" ]] || [[ $stage == "Integration_DayPatient" ]] || [[ $stage == "Clustree_Integration_DayPatient" ]] || [[ $stage == "QCfigures_Integration_DayPatient" ]] || [[ $stage == "Integration_Day" ]] || [[ $stage == "Integration_Site" ]] || [[ $stage == "Integration_Site_Regress_MtRb" ]] || [[ $stage == "Integration_Day_then_Individual" ]]  || [[ $stage == "Clustree_Integration_Day_then_Individual" ]] || [[ $stage == "QCfigures_Integration_Day_then_Individual" ]]
then
    PE=1
    T=1-1
    MEM=100G
elif [[ $stage == "Multi_Resolution_no_covatiates" ]] || [[ $stage == "Multi_Resolution_PoolRegressed" ]] || [[ $stage == "Multi_Resolution_PoolDayRegressed" ]] || [[ $stage == "Multi_Resolution_Harmony_Pool" ]] || [[ $stage == "Multi_Resolution_Harmony_Day" ]] || [[ $stage == "Multi_Resolution_Harmony_DayPatient" ]] || [[ $stage == "Multi_Resolution_Integration_DayPatient" ]] || [[ $stage == "Multi_Resolution_Integration_Day_then_Individual" ]]
then
    T=1-25
    # T=3,7
    PE=2
    MEM=20G
elif [[ $stage == "SCTnormalization_SiteMtRbRegressed" ]] 
then
    PE=2
    T=1-1
    MEM=100G
fi

qsub -S /bin/bash \
        -q short.q \
        -r yes \
        -l mem_requested=$MEM \
        -l tmp_requested=$MEM \
        -N Seurat_$stage \
        -pe smp $PE \
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