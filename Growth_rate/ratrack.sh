#!/bin/bash

RATRACK=/directflow/SCCGGroupShare/projects/DrewNeavin/software/ratrack/
SNAKEFILE=/directflow/SCCGGroupShare/projects/DrewNeavin/software/ratrack/Snakefile


cd $RATRACK
mkdir logs

conda activate ratrack

# ### First test with demo ###
# snakemake results/minimal.pdf results/minimal.fit.csv --dryrun
# snakemake results/kcl22.pdf results/kcl22.fit.csv


# snakemake results/TOB0421_ratrack.pdf results/TOB0421_ratrack.fit.csv --dryrun
# nohup snakemake results/TOB0421_ratrack.pdf results/TOB0421_ratrack.fit.csv --cores 5 > logs/nohup_`date +%Y-%m-%d.%H:%M:%S`_TOB0421.log &

# snakemake results/FSA0006_ratrack.pdf results/FSA0006_ratrack.fit.csv --dryrun
# nohup snakemake results/FSA0006_ratrack.pdf results/FSA0006_ratrack.fit.csv --cores 5 > logs/nohup_`date +%Y-%m-%d.%H:%M:%S`_FSA0006.log &

# snakemake results/MBE1006_ratrack.pdf results/MBE1006_ratrack.fit.csv --dryrun
# nohup snakemake results/MBE1006_ratrack.pdf results/MBE1006_ratrack.fit.csv --cores 5 > logs/nohup_`date +%Y-%m-%d.%H:%M:%S`_MBE1006.log &


### Try with original counts + dillution info ###
snakemake results/TOB0421_original_counts.pdf results/TOB0421_original_counts.fit.csv --dryrun
nohup snakemake results/TOB0421_original_counts.pdf results/TOB0421_original_counts.fit.csv --cores 5 > logs/nohup_`date +%Y-%m-%d.%H:%M:%S`_TOB0421.log &

snakemake results/FSA0006_original_counts.pdf results/FSA0006_original_counts.fit.csv --dryrun
nohup snakemake results/FSA0006_original_counts.pdf results/FSA0006_original_counts.fit.csv --cores 5 > logs/nohup_`date +%Y-%m-%d.%H:%M:%S`_FSA0006.log &

snakemake results/MBE1006_original_counts.pdf results/MBE1006_original_counts.fit.csv --dryrun
nohup snakemake results/MBE1006_original_counts.pdf results/MBE1006_original_counts.fit.csv --cores 5 > logs/nohup_`date +%Y-%m-%d.%H:%M:%S`_MBE1006.log &





snakemake results/TOB0421_ratrack.pdf results/TOB0421_ratrack.fit.csv --dryrun --unlock
