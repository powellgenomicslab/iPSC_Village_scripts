import pandas as pd


genes_file = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/multi-passage/deboever_imputed_overlapping_filtered_header_pruned_snp_gene.tsv"
# genes_file_kilpinen = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/KilpinenOverlap/kilpinen_imputed_overlapping_filtered_header_pruned_snp_gene.tsv"
genes = pd.read_csv(genes_file, sep = "\t")
# genes_kilpenen = pd.read_csv(genes_file_kilpinen, sep = "\t")
# genes_kilpenen = genes_kilpenen.iloc[1]

rule all:
    input:
        expand("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/deboever/gene_separate/beds/{gene}_deboever_eQTL_results.bed", gene = genes.gene_id),
        # expand("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/kilpinen/gene_separate/beds/{gene}_kilpinen_eQTL_results.bed", gene = genes_kilpenen.gene_id),


rule eqtl:
    input:
        bed = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/multi-passage/deboever_imputed_overlapping_filtered_header_pruned.bed"
    output:
        "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/deboever/gene_separate/beds/{gene}_deboever_eQTL_results.bed"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 4,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 4
    threads: 4
    params:
        script = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/eQTL_check/multi-passage/test_eQTL.R",
        outdir="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/deboever/gene_separate/",
        snp=lambda wildcards: genes.ID[genes.gene_id == wildcards.gene],
    log:
    shell:
        """        
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/baseR402/bin/Rscript {params.script} {wildcards.gene} {params.outdir} {input.bed}
        """


rule eqtl_kilpinen:
    input:
        bed = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/KilpinenOverlap/kilpinen_imputed_overlapping_filtered_header_pruned.bed"
    output:
        "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/kilpinen/gene_separate/beds/{gene}_kilpinen_eQTL_results.bed"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 4,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 4
    threads: 4
    params:
        script = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/eQTL_check/test_eQTL_kilpinen.R",
        outdir="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/kilpinen/gene_separate/",
        snp=lambda wildcards: genes.ID_ref_alt[genes.gene_id == wildcards.gene],
    log:
    shell:
        """        
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/baseR402/bin/Rscript {params.script} {wildcards.gene} {params.outdir} {input.bed}
        """
