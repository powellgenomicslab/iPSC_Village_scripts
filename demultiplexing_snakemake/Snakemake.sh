#!/bin/bash

export SINGULARITY_BINDPATH="/directflow/"
SNAKEFILE="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/demultiplexing_snakemake/Snakefile"
LOGS="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/demultiplexing_snakemake/logs"
mkdir -p $LOGS

nohup snakemake --snakefile $SNAKEFILE \
    --rerun-incomplete --jobs 100 --use-singularity --restart-times 4 --keep-going \
    --cluster "qsub -S /bin/bash -q short.q -r yes -pe smp {threads} -l tmp_requested={resources.disk_per_thread_gb}G -l mem_requested={resources.mem_per_thread_gb}G -e $LOGS -o $LOGS -j y -V -m e -M d.neavin@garvan.org.au -l hostname='zeta*|beefy*'" \
    > $LOGS/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &

snakemake --snakefile $SNAKEFILE \
    --rerun-incomplete --jobs 100 --use-singularity --restart-times 4 --dryrun --reason --unlock