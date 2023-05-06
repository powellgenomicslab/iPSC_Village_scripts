

SNAKEFILE="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/eQTL_check/uni_village/village/test_eQTL.smk"
LOG="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/logs/"

mkdir -p $LOG
cd $LOG

snakemake \
  --snakefile $SNAKEFILE \
  --dryrun \
  --cores 1 \
  --quiet \
  --unlock


snakemake \
  --snakefile $SNAKEFILE \
  --dryrun \
  --cores 1 \
  --reason > jobs2run.txt



nohup \
  snakemake \
    --snakefile $SNAKEFILE \
    --jobs 200 \
    --use-singularity \
    --restart-times 1 \
    --keep-going \
    --cluster \
        "qsub -S /bin/bash \
        -q short.q \
        -r yes \
        -pe smp {threads} \
        -l tmp_requested={resources.disk_per_thread_gb}G \
        -l mem_requested={resources.mem_per_thread_gb}G \
        -e $LOG \
        -o $LOG \
        -j y \
        -V" \
  > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &


snakemake \
  --snakefile $SNAKEFILE \
  --dryrun \
  --cores 1 \
  --quiet \
  --unlock


snakemake \
  --snakefile $SNAKEFILE \
  --dryrun \
  --cores 1 \
  --quiet \
  --cleanup-metadata \
  /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review/gene_separated/icc/ENSG00000116001_icc.rds \
  /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review/gene_separated/fit_models/ENSG00000116001_fitted_models.rds \
  /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review/gene_separated/icc/ENSG00000084112_icc.rds \
  /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review/gene_separated/fit_models/ENSG00000084112_fitted_models.rds \
  /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review/gene_separated/icc/ENSG00000130706_icc.rds \
  /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review/gene_separated/fit_models/ENSG00000130706_fitted_models.rds






ENSG00000140481
ENSG00000158552
ENSG00000167004
ENSG00000184281
ENSG00000185437
ENSG00000185875
ENSG00000196109
ENSG00000197536
ENSG00000215910
ENSG00000215883
ENSG00000258297
ENSG00000259488
ENSG00000274877


ENSG00000158555
ENSG00000167005
ENSG00000184292
ENSG00000184281
ENSG00000185437
ENSG00000185875
ENSG00000196109
ENSG00000197536
ENSG00000215866
ENSG00000258289
ENSG00000259485
ENSG00000274828


rm ENSG00000167011_fitted_models.rds
rm ENSG00000184304_fitted_models.rds
rm ENSG00000185477_fitted_models.rds
rm ENSG00000185900_fitted_models.rds
rm ENSG00000196132_fitted_models.rds
rm ENSG00000197566_fitted_models.rds
rm ENSG00000217930_fitted_models.rds
rm ENSG00000258593_fitted_models.rds
rm ENSG00000259715_fitted_models.rds
rm ENSG00000275074_fitted_models.rds
