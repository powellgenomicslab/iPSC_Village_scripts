
SNAKEFILE="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/post_review_all_cells/variance_partitioning_all_cells.smk"
LOG="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partitioning_all_cells/gene_separated/logs/"

mkdir -p $LOG


conda activate wg1_snakemake


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



rm ENSG00000047617_icc.rds
rm ENSG00000089250_icc.rds
rm ENSG00000124602_icc.rds
rm ENSG00000130052_icc.rds
rm ENSG00000136099_icc.rds
rm ENSG00000143412_icc.rds
rm ENSG00000147124_icc.rds
rm ENSG00000158104_icc.rds
rm ENSG00000160188_icc.rds
rm ENSG00000163393_icc.rds
rm ENSG00000163947_icc.rds
rm ENSG00000167083_icc.rds
rm ENSG00000170113_icc.rds
rm ENSG00000170961_icc.rds
rm ENSG00000172748_icc.rds
rm ENSG00000184916_icc.rds
rm ENSG00000185386_icc.rds
rm ENSG00000186141_icc.rds
rm ENSG00000188536_icc.rds
rm ENSG00000189164_icc.rds
rm ENSG00000230606_icc.rds
rm ENSG00000242687_icc.rds
rm ENSG00000243444_icc.rds
rm ENSG00000253438_icc.rds
rm ENSG00000253731_icc.rds
rm ENSG00000258813_icc.rds
rm ENSG00000258944_icc.rds
rm ENSG00000260528_icc.rds
rm ENSG00000272030_icc.rds
rm ENSG00000272973_icc.rds
rm ENSG00000274367_icc.rds
rm ENSG00000275580_icc.rds




rm ENSG00000047617_fitted_models.rds
rm ENSG00000089250_fitted_models.rds
rm ENSG00000124602_fitted_models.rds
rm ENSG00000130052_fitted_models.rds
rm ENSG00000136099_fitted_models.rds
rm ENSG00000143412_fitted_models.rds
rm ENSG00000147124_fitted_models.rds
rm ENSG00000158104_fitted_models.rds
rm ENSG00000160188_fitted_models.rds
rm ENSG00000163393_fitted_models.rds
rm ENSG00000163947_fitted_models.rds
rm ENSG00000167083_fitted_models.rds
rm ENSG00000170113_fitted_models.rds
rm ENSG00000170961_fitted_models.rds
rm ENSG00000172748_fitted_models.rds
rm ENSG00000184916_fitted_models.rds
rm ENSG00000185386_fitted_models.rds
rm ENSG00000186141_fitted_models.rds
rm ENSG00000188536_fitted_models.rds
rm ENSG00000189164_fitted_models.rds
rm ENSG00000230606_fitted_models.rds
rm ENSG00000242687_fitted_models.rds
rm ENSG00000243444_fitted_models.rds
rm ENSG00000253438_fitted_models.rds
rm ENSG00000253731_fitted_models.rds
rm ENSG00000258813_fitted_models.rds
rm ENSG00000258944_fitted_models.rds
rm ENSG00000260528_fitted_models.rds
rm ENSG00000272030_fitted_models.rds
rm ENSG00000272973_fitted_models.rds
rm ENSG00000274367_fitted_models.rds
rm ENSG00000275580_fitted_models.rds
