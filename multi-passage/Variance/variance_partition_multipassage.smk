import pandas as pd


genes_file = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/preQC/genes_1pct_expressing.tsv"
genes = pd.read_csv(genes_file, sep = "\t")
# genes = genes.iloc[1]


rule all:
    input:
        # expand("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance/gene_separated/fit_models/{gene}_fitted_models.rds", gene = genes.Gene),
        # expand("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_ncov/gene_separated/fit_models/{gene}_fitted_models.rds", gene = genes.Gene),
        expand("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_integratedSCT/gene_separated/fit_models/{gene}_fitted_models.rds", gene = genes.Gene),
        expand("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_integratedSCT/gene_separated/icc/{gene}_icc.rds", gene = genes.Gene),


rule partition_variance:
    input:
       "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/preQC/seurat_joint_SCT_1pct_expressing.rds"
    output:
        "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance/gene_separated/icc/{gene}_icc.rds",
        "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance/gene_separated/fit_models/{gene}_fitted_models.rds"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 16,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 16
    threads: 4
    params:
        script = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/multi-passage/Variance/variance_partition_multipassage.R",
        out_icc="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance/gene_separated/icc/",
        out_model="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance/gene_separated/fit_models/",
        out_resids="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance/gene_separated/residuals4qtl/",
        out_icc_interaction = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance/gene_separated/icc_interaction/",
        out_model_interaction = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance/gene_separated/icc_interaction/"
    log:
    shell:
        """
        mkdir -p {params.out_icc_interaction}
        mkdir -p {params.out_model_interaction}
        mkdir -p {params.out_resids}
        
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/baseR402/bin/Rscript {params.script} {params.out_icc_interaction} {params.out_icc} {params.out_model_interaction} {params.out_model} {params.out_resids} {wildcards.gene}
        """


rule partition_variance_ncov:
    input:
        seurat = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/preQC/seurat_joint_SCT_1pct_expressing.rds"
    output:
        "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_ncov/gene_separated/icc/{gene}_icc.rds",
        "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_ncov/gene_separated/icc_interaction/{gene}_icc.rds",
        "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_ncov/gene_separated/fit_models/{gene}_fitted_models.rds"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 16,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 16
    threads: 4
    params:
        script = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/multi-passage/Variance/variance_partition_multipassage_integratedSCT_Ncov.R",
        out_icc="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_ncov/gene_separated/icc/",
        out_model="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_ncov/gene_separated/fit_models/",
        out_resids="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_ncov/gene_separated/residuals4qtl/",
        out_icc_interaction = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_ncov/gene_separated/icc_interaction/",
        out_model_interaction = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_ncov/gene_separated/icc_interaction/"
    log:
    shell:
        """
        mkdir -p {params.out_icc_interaction}
        mkdir -p {params.out_model_interaction}
        mkdir -p {params.out_resids}
        
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/baseR402/bin/Rscript {params.script} {params.out_icc_interaction} {params.out_icc} {params.out_model_interaction} {params.out_model} {params.out_resids} {wildcards.gene}
        """
        

### Tried manually with just CHCHD2 and the Ncov impacts the ability to get the line effect so cannot run this way (could downsample some of the pools to be more even sizes)
rule partition_variance_integratedSCT:
    input:
        seurat =  "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/preQC/time-integrated_filtered_seurat_1pct_expressing.rds"
    output:
        "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_integratedSCT/gene_separated/icc/{gene}_icc.rds",
        "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_integratedSCT/gene_separated/fit_models/{gene}_fitted_models.rds"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 8,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 8
    threads: 4
    params:
        script = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/multi-passage/Variance/variance_partition_multipassage_integratedSCT.R",
        out_icc="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_integratedSCT/gene_separated/icc/",
        out_model="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_integratedSCT/gene_separated/fit_models/",
        out_resids="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_integratedSCT/gene_separated/residuals4qtl/",
        out_icc_interaction = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_integratedSCT/gene_separated/icc_interaction/",
        out_model_interaction = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_integratedSCT/gene_separated/icc_interaction/"
    log:
    shell:
        """
        mkdir -p {params.out_icc_interaction}
        mkdir -p {params.out_model_interaction}
        mkdir -p {params.out_resids}
        
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/baseR402/bin/Rscript {params.script} {params.out_icc_interaction} {params.out_icc} {params.out_model_interaction} {params.out_model} {params.out_resids} {wildcards.gene}
        """