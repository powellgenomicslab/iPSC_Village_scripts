import pandas as pd


genes_file = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/uni_village/cryopreserved/village/deboever_imputed_overlapping_filtered_header_pruned_snp_gene.tsv"
genes = pd.read_csv(genes_file, sep = "\t")
# genes = genes.iloc[1]

rule all:
    input:
        expand("/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/uni_village/cryopreserved/village/gene_separate/beds/{gene}_deboever_eQTL_results.bed", gene = genes.gene_id),


rule eqtl:
    input:
        bed = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/uni_village/cryopreserved/uniculture/deboever_imputed_overlapping_filtered_header_pruned.bed"
    output:
        "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/uni_village/cryopreserved/village/gene_separate/beds/{gene}_deboever_eQTL_results.bed"
    resources:
        mem_per_thread_gb = lambda wildcards, attempt: attempt * 4,
        disk_per_thread_gb = lambda wildcards, attempt: attempt * 4
    threads: 4
    params:
        script = "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/eQTL_check/uni_village/cryopreserved/village/test_eQTL.R",
        outdir="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/uni_village/cryopreserved/village/gene_separate/",
    log:
    shell:
        """        
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/baseR402/bin/Rscript {params.script} {wildcards.gene} {params.outdir} {input.bed}
        """

