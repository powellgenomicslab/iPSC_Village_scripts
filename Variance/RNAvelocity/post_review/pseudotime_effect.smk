import pandas as pd


genes_file = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/post_review/data/seurat_integrated_all_times_clustered_1pct_expressing_pseudotime.tsv"
genes = pd.read_csv(genes_file, sep = "\t")
# genes = genes.iloc[1:1000]


rule all:
    input:
        expand( "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/post_review/results/gene_separated/icc/{gene}_icc.rds", gene = genes.Gene)


rule partition_variance:
    input:
        seurat = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/post_review/data/seurat_integrated_all_times_clustered_1pct_expressing_pseudotime.rds"
    output:
        "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/post_review/results/gene_separated/icc/{gene}_icc.rds"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 64,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 64
    threads: 1
    params:
        script = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/Variance/RNAvelocity/post_review/pseudotime_effect.R",
        out_icc="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/post_review/results/gene_separated/icc/",
        out_icc_interaction = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/post_review/results/gene_separated/icc_interaction/",
        out_model_interaction = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/post_review/results/gene_separated/icc_interaction/",
        out_plot = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/post_review/results/gene_separated/plots/",
        out_effects = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/Variance/RNAvelocity/post_review/results/gene_separated/effect_betas/"
    log:
    shell:
        """
        mkdir -p {params.out_icc_interaction}
        mkdir -p {params.out_model_interaction}
        mkdir -p {params.out_plot}
        mkdir -p {params.out_effects}
        
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/baseR402/bin/Rscript {params.script} {params.out_icc_interaction} {params.out_icc} {params.out_plot} {wildcards.gene} {params.out_effects}
        """


