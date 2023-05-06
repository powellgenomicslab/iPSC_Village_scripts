#!/bin/bash

RATRACK=/directflow/SCCGGroupShare/projects/DrewNeavin/software/ratrack/
SNAKEFILE=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Growth_rate/Snakefile_edited_more_times.smk

OUTDIR=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/growth_rate/18_line_village/ratrack/data/more_growth_rate_points/
LOGS=$OUTDIR/logs


mkdir $LOGS

conda activate ratrack

cd /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Growth_rate


### Try with original counts + dillution info ###
snakemake --snakefile $SNAKEFILE \
            --dryrun \
            --cores 1 \
            --use-conda \
            --rerun-incomplete

# snakemake --snakefile $SNAKEFILE \
#             --cores 5 \
#             --use-conda \
#             --keep-going

snakemake \
    --snakefile $SNAKEFILE \
    --rerun-incomplete \
    --jobs 20 \
    --use-conda \
    --keep-going \
    --cluster \
        "qsub -S /bin/bash \
        -q short.q \
        -r yes \
        -pe smp 4 \
        -l tmp_requested=8G \
        -l mem_requested=8G \
        -e $LOGS \
        -o $LOGS \
        -j y \
        -V"



        