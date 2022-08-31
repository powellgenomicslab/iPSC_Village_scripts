#!/bin/bash

OUTDIR="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/multi-passage"
mkdir -p $OUTDIR

kilpinen_genes="$OUTDIR/../KilpinenOverlap/gene_snp_list.tsv"
deboever_genes="$OUTDIR/../KilpinenOverlap/deboever_gene_snp_list.tsv"
variance_genes="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/Variance/variance_integratedSCT/gene_separated/residuals4qtl"


ls $variance_genes | sed "s/_residuals4qtl.rds//g" | sort -u > $OUTDIR/residual_genes.tsv
awk '{print $2}' $kilpinen_genes | sort -u > $OUTDIR/kilpinen_genes.tsv
awk '{print $2}' $deboever_genes | sort -u > $OUTDIR/deboever_genes.tsv


comm -12 $OUTDIR/residual_genes.tsv head $OUTDIR/kilpinen_genes.tsv > $OUTDIR/genes_residual_and_kilpinen.tsv ## NONE
comm -12 $OUTDIR/residual_genes.tsv $OUTDIR/deboever_genes.tsv > $OUTDIR/genes_residual_and_deboever.tsv


### Get finalized list of genes that were eQTLs and demonstrate line effects
grep -F -f $OUTDIR/genes_residual_and_deboever.tsv $deboever_genes > $OUTDIR/finalized_deboever_gene_snp_list.tsv



### Filter vcf for just these snps ###
awk '{print $1}' $OUTDIR/finalized_deboever_gene_snp_list.tsv | sort -u > $OUTDIR/finalized_deboever_snps.tsv

grep "#" $OUTDIR/nona_cardiac_multiome_Filtered_INFO_0.4_MAF0.05_complete_cases_snps_filtered_diff_genotypes.vcf >$OUTDIR/deboever_finalized_snps.vcf
grep -f $OUTDIR/finalized_deboever_snps.tsv $OUTDIR/nona_cardiac_multiome_Filtered_INFO_0.4_MAF0.05_complete_cases_snps_filtered_diff_genotypes.vcf >> $OUTDIR/deboever_finalized_snps.vcf

