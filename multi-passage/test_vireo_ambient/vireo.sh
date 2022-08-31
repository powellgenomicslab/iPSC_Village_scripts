#!/bin/bash




singularity exec --bind /directflow $SIF vireo -c $VIREO_INDIR -d $VIREO_INDIR/donor_subset.vcf -o $VIREO_OUTDIR -t $FIELD --callAmbientRNAs