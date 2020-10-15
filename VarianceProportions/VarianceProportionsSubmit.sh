#!/bin/bash

# stage="SCTnormalizedReplicate"
# stage="SCTnormalized"
# stage="scaterVariance"
# stage="scaterVarianceALL"
# stage="scaterVariance_RbMtSiteCov"
# stage="scaterVariance_SeparatedSite"
# stage="scaterVariance_SeparatedSite_MtRb_regressed"
stage="scaterVariance_SeparatedReplicate"
# stage="scaterVariance_SeparatedSiteReplicate"
# stage="scaterVariance_SeparatedSiteReplicate_RbMtRegressed"


DIR="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/"
OUT="$DIR/output/VarianceProportions/$stage/"
LOG="$OUT/logs/"
PIPELINE="$DIR/scripts/VarianceProportions/VarianceProportions.sh"
ARGS="$OUT/ARGS/"

mkdir -p $LOG
mkdir -p $OUT
mkdir -p $ARGS

export PATH=$PATH

today=`date +%Y-%m-%d.%H:%M:%S`
ARGfile="$ARGS/arguments$stage$today.txt"

echo $stage > $ARGfile
echo $OUT >> $ARGfile
echo $DIR >> $ARGfile

cd $OUT
if [[ $stage == "scaterVariance_SeparatedSite" ]] || [[ $stage == "scaterVariance_SeparatedSite_MtRb_regressed" ]]
then
        PE=5
        MEM=20G
        T=1-3
elif [[ $stage == "scaterVariance_SeparatedReplicate" ]]
then
        PE=5
        MEM=5G
        T=6,10
elif [[ $stage == "scaterVariance_SeparatedSiteReplicate" ]] || [[ $stage == "scaterVariance_SeparatedSiteReplicate_RbMtRegressed" ]]
then
        PE=5
        MEM=5G
        T=1-9
else
        PE=5
        MEM=100G
        T=1-1
fi

qsub -S /bin/bash \
        -q short.q \
        -r yes \
        -l mem_requested=$MEM \
        -l tmp_requested=$MEM \
        -N var_prop_$stage \
        -pe smp $PE \
        -cwd \
        -m e \
        -hold_jid 360112 \
        -M d.neavin@garvan.org.au \
        -j y \
        -e $LOG \
        -o $LOG \
        -V \
        -t $T \
        -v DIR=$DIR,ARGfile=$ARGfile,stage=$stage,today=$today,OUT=$OUT,LOG=$LOG,PIPELINE=$PIPELINE,ARGS=$ARGS \
        -C '' $PIPELINE

