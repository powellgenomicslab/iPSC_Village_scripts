#!/bin/bash

SIF=/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/workflow_testing/v2.0.0/Demuxafy.sif

DIR_ONEK=/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/output/Consortium/benchmark_rerun
DIR_FIBRO=/directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/output/fibroblasts
DIR_VILLAGE=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/CombinedResults
DIR_POAG=/directflow/SCCGGroupShare/projects/annsen/analysis/POAG_scRNA/demuxlet

OUTDIR=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Pooling_variance

mkdir -p $OUTDIR


for directory in `ls -d $DIR_ONEK/OneK1K_scRNA_Sample* | grep -v "_V1"`
do
    pool=`basename $directory`
    mkdir -p $OUTDIR/$pool
    singularity exec --bind /directflow $SIF Combine_Results.R \
        -o $OUTDIR/$pool/combined_assignments \
        -d $DIR_ONEK/$pool/popscle/demuxlet/demuxletOUT_impute_vars.best \
        -f $DIR_ONEK/$pool/popscle/freemuxlet/freemuxletOUT.clust1.samples.gz \
        -u $DIR_ONEK/$pool/souporcell/clusters.tsv \
        -v $DIR_ONEK/$pool/vireo/results/donor_ids.tsv \
        --method "MajoritySinglet"
done


for directory in `ls -d $DIR_FIBRO/scFibroblast_EQTL_Sample*`
do
    pool=`basename $directory`
    mkdir -p $OUTDIR/$pool
    singularity exec --bind /directflow $SIF Combine_Results.R \
        -o $OUTDIR/$pool/combined_assignments \
        -d $DIR_FIBRO/$pool/popscle/demuxlet/demuxletOUT_impute_vars.best \
        -f $DIR_FIBRO/$pool/popscle/freemuxlet/freemuxletOUT.clust1.samples.gz \
        -u $DIR_FIBRO/$pool/souporcell/clusters.tsv \
        -v $DIR_FIBRO/$pool/vireo/results/donor_ids.tsv \
        --method "MajoritySinglet"
done




for directory in `ls -d $DIR_VILLAGE/*`
do
    pool=`basename $directory`
    mkdir -p $OUTDIR/$pool
    singularity exec --bind /directflow $SIF Combine_Results.R \
        -o $OUTDIR/$pool/combined_assignments \
        -d /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/popscle/demuxlet/$pool/Imputed_GPdemuxletOUT.best \
        -f /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/popscle/freemuxlet/$pool/freemuxletOUT.clust1.samples.gz \
        -s /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scSplit/$pool/scSplit_result.csv \
        -t /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/DoubletDetection/$pool/DoubletDetection_results.txt \
        -u /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/souporcell/$pool/clusters.tsv \
        -v /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/vireo/$pool/results/donor_ids.tsv \
        --method "MajoritySinglet"
done



