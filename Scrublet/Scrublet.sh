#!/bin/bash

SAMPLE=`head -n $SGE_TASK_ID $SAMPLEINFO | tail -n 1 | awk '{print $1}'`
INPUT=`head -n $SGE_TASK_ID $SAMPLEINFO | tail -n 1 | awk '{print $2}'`

INPUTDIR="$INPUT/Expression_200128_A00152_0196_BH3HNFDSXY/GE/$SAMPLE//outs/filtered_feature_bc_matrix/"

OUTPUT=$OUT/$SAMPLE

mkdir -p $OUTPUT

source activate scrublet
python $DIR/scripts/Scrublet/Scrublet.py $INPUTDIR $OUTPUT $VAR $THRESHOLD
source deactivate