import pandas as pd


genes_file = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review_cryo/seurat_integrated_Sydney_1pct_expressing_genes.tsv"
genes = pd.read_csv(genes_file, sep = "\t")


rule all:
    input:
        expand("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review_cryo/gene_separated/fit_models/{gene}_fitted_models.rds", gene = genes.Gene)


rule partition_variance:
    input:
        seurat = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Distribution_tests/seurat_integrated_all_times_clustered_1pct_expressing.rds"
    output:
        "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review_cryo/gene_separated/icc/{gene}_icc.rds",
        "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review_cryo/gene_separated/fit_models/{gene}_fitted_models.rds"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 16,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 16
    threads: 4
    params:
        script = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/post_review_cryo/variance_partition_post_review_cryo.R",
        out_icc="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review_cryo/gene_separated/icc/",
        out_model="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review_cryo/gene_separated/fit_models/",
        out_resids="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review_cryo/gene_separated/residuals4qtl/"
    log:
    shell:
        """
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/baseR402/bin/Rscript {params.script} {params.out_icc} {params.out_model} {params.out_resids} {wildcards.gene}
        """
        