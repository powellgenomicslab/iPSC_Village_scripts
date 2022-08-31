import pandas as pd


genes_file = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/KilpinenOverlap/deboever_imputed_overlapping_filtered_header_pruned_snp_gene.tsv"
genes_file_kilpinen = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/KilpinenOverlap/kilpinen_imputed_overlapping_filtered_header_pruned_snp_gene.tsv"
genes = pd.read_csv(genes_file, sep = "\t")
# genes = genes.iloc[1]
genes_kilpenen = pd.read_csv(genes_file_kilpinen, sep = "\t")
# genes_kilpenen = genes_kilpenen.iloc[1]

rule all:
    input:
        expand("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/eQTL_village_interaction/deboever/results/{gene}_snpXvillage_interactions.tsv", gene = genes.gene_id),
        expand("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/eQTL_village_interaction/kilpinen/results/{gene}_snpXvillage_interactions.tsv", gene = genes_kilpenen.gene_id),


rule interaction_deboever:
    input:
        bed = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/KilpinenOverlap/deboever_imputed_overlapping_filtered_header_pruned.bed"
    output:
        "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/eQTL_village_interaction/deboever/results/{gene}_snpXvillage_interactions.tsv",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 8,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 8
    threads: 4
    params:
        script = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/eQTL_check/eQTL_village_interaction_deboever.R",
        outdir="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/eQTL_village_interaction/deboever/results/",
        snp=lambda wildcards: genes.ID_ref_alt[genes.gene_id == wildcards.gene],
        datadir = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/eQTL_village_interaction/deboever/data/"
    log:
    shell:
        """        
        mkdir -p {params.datadir}
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/baseR402/bin/Rscript {params.script} {wildcards.gene} {params.outdir} {input.bed} {params.datadir}
        """



rule interaction_kilpinen:
    input:
        bed = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/KilpinenOverlap/kilpinen_imputed_overlapping_filtered_header_pruned.bed"
    output:
        "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/eQTL_village_interaction/kilpinen/results/{gene}_snpXvillage_interactions.tsv",
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 8,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 8
    threads: 4
    params:
        script = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/eQTL_check/eQTL_village_interaction_kilpinen.R",
        outdir="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/eQTL_village_interaction/kilpinen/results/",
        snp=lambda wildcards: genes.ID_ref_alt[genes.gene_id == wildcards.gene],
        datadir = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/eQTL_village_interaction/kilpinen/data/"
    log:
    shell:
        """  
        mkdir -p {params.datadir}
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/baseR402/bin/Rscript {params.script} {wildcards.gene} {params.outdir} {input.bed} {params.datadir}
        """
