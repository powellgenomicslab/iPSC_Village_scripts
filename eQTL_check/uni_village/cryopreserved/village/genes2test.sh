#!/bin/bash

DIR="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/uni_village"
OUTDIR="$DIR/cryopreserved/village"

mkdir -p $OUTDIR

deboever_genes="$DIR/deboever_gene_snp_list.tsv"
variance_genes="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review_cryo2/village/gene_separated/residuals4qtl"


ls $variance_genes | sed "s/_residuals4qtl.rds//g" | sort -u > $OUTDIR/residual_genes.tsv
awk '{print $2}' $deboever_genes | sort -u > $OUTDIR/deboever_genes.tsv


comm -12 $OUTDIR/residual_genes.tsv $OUTDIR/deboever_genes.tsv > $OUTDIR/genes_residual_and_deboever.tsv


### Get finalized list of genes that were eQTLs and demonstrate line effects
grep -F -f $OUTDIR/genes_residual_and_deboever.tsv $deboever_genes > $OUTDIR/finalized_deboever_gene_snp_list.tsv



### Filter vcf for just these snps ###
awk '{print $1}' $OUTDIR/finalized_deboever_gene_snp_list.tsv | sort -u > $OUTDIR/finalized_deboever_snps.tsv


conda activate vcftools

    vcftools --vcf $DIR/merged_imputed_AllChrs_iPSC_R2_0.3_filtered_diff_genotypes.vcf --snps $OUTDIR/finalized_deboever_snps.tsv --recode --recode-INFO-all --out $OUTDIR/deboever_finalized_snps

conda deactivate


# ### Get 1000G ref and alt allele annotations
# vcf1000g="/directflow/SCCGGroupShare/projects/DrewNeavin/References/GRCh37SNPvcfs1000genomes/MergedAutosomesFilteredGenes.recode.MAF0.01.vcf"

# conda activate bcftools

#     bgzip -c $OUTDIR/deboever_finalized_snps.recode.vcf > $OUTDIR/deboever_finalized_snps.recode.vcf.gz
#     tabix -p vcf $OUTDIR/deboever_finalized_snps.recode.vcf.gz

#     bgzip -c $vcf1000g > $OUTDIR/MergedAutosomesFilteredGenes_MAF0.01.1000g_hg19.vcf.gz
#     tabix -p vcf $OUTDIR/MergedAutosomesFilteredGenes_MAF0.01.1000g_hg19.vcf.gz

#     bcftools isec $OUTDIR/deboever_finalized_snps.recode.vcf.gz $OUTDIR/MergedAutosomesFilteredGenes_MAF0.01.1000g_hg19.vcf.gz -o $OUTDIR/overlap_1000g.vcf.gz -O vz -p $OUTDIR

# conda deactivate


