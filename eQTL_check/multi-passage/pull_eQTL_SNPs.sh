#!/bin/bash


vcf=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/multi-passage/nona_cardiac_multiome_Filtered_INFO_0.4_MAF0.05_complete_cases.recode.vcf
kilpinen=/directflow/SCCGGroupShare/projects/DrewNeavin/References/iPSC_eQTLs/KilpineniPSCeQTLs.txt
deboever=/directflow/SCCGGroupShare/projects/DrewNeavin/References/iPSC_eQTLs/DeBoeveriPSCeQTLs.txt

OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/References/iPSC_eQTLs
OUT_TMP=/directflow/SCCGGroupShare/projects/DrewNeavin/References/iPSC_eQTLs/tmp

mkdir -p $OUT_TMP

INTERSECT_OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/multi-passage

mkdir -p $INTERSECT_OUT



### Filter on individuals and MAF for our lines ###
sed '/<CN/d' $vcf > $INTERSECT_OUT/nona_cardiac_multiome_Filtered_INFO_0.4_MAF0.05_complete_cases_snps.vcf
sed -i '/<INV>/d' $INTERSECT_OUT/nona_cardiac_multiome_Filtered_INFO_0.4_MAF0.05_complete_cases_snps.vcf



conda activate vcftools

    vcftools --vcf $INTERSECT_OUT/nona_cardiac_multiome_Filtered_INFO_0.4_MAF0.05_complete_cases_snps.vcf --recode --mac 1 --out $INTERSECT_OUT/nona_cardiac_multiome_Filtered_INFO_0.4_MAF0.05_complete_cases_snps_filtered

    bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES || COUNT(GT="RA")=N_SAMPLES' $INTERSECT_OUT/nona_cardiac_multiome_Filtered_INFO_0.4_MAF0.05_complete_cases_snps_filtered.recode.vcf > $INTERSECT_OUT/nona_cardiac_multiome_Filtered_INFO_0.4_MAF0.05_complete_cases_snps_filtered_diff_genotypes.vcf

conda deactivate

### Redo intersection with these SNPs ##$#
conda activate bedtools

    bedtools intersect -a $OUT/kilpinen.bed -b $INTERSECT_OUT/nona_cardiac_multiome_Filtered_INFO_0.4_MAF0.05_complete_cases_snps_filtered_diff_genotypes.vcf -wa -wb > $INTERSECT_OUT/kilpinen_imputed_overlapping_filtered.bed

conda deactivate



head -n 1 $OUT/kilpinen.bed > $INTERSECT_OUT/kilpinen_header_bed.tsv
grep "#CHROM" $INTERSECT_OUT/nona_cardiac_multiome_Filtered_INFO_0.4_MAF0.05_complete_cases_snps_filtered_diff_genotypes.vcf > $INTERSECT_OUT/kilpinen_header_vcf.tsv
paste -d"\t" $INTERSECT_OUT/kilpinen_header_bed.tsv $INTERSECT_OUT/kilpinen_header_vcf.tsv > $INTERSECT_OUT/kilpinen_combined_header.tsv

cat $INTERSECT_OUT/kilpinen_combined_header.tsv > $INTERSECT_OUT/kilpinen_imputed_overlapping_filtered_header.bed
cat $INTERSECT_OUT/kilpinen_imputed_overlapping_filtered.bed >> $INTERSECT_OUT/kilpinen_imputed_overlapping_filtered_header.bed





### Use this file to test for trends in gene
awk 'BEGIN{FS=OFS="\t"}{print($24, $7)}' $INTERSECT_OUT/kilpinen_imputed_overlapping_filtered_header.bed | sed 's/\..*//g' > $INTERSECT_OUT/gene_snp_list.tsv

sed -i '1i snp\tgene' $INTERSECT_OUT/gene_snp_list.tsv







##### DEBOEVER DATA #####
### Intersect with these SNPs ##$#
conda activate bedtools

    bedtools intersect -a $OUT/deboever.bed -b $INTERSECT_OUT/nona_cardiac_multiome_Filtered_INFO_0.4_MAF0.05_complete_cases_snps_filtered_diff_genotypes.vcf -wa -wb > $INTERSECT_OUT/deboever_imputed_overlapping_filtered.bed

conda deactivate


head -n 1 $OUT/deboever.bed > $INTERSECT_OUT/deboever_header_bed.tsv
grep "#CHROM" $INTERSECT_OUT/nona_cardiac_multiome_Filtered_INFO_0.4_MAF0.05_complete_cases_snps_filtered_diff_genotypes.vcf > $INTERSECT_OUT/deboever_header_vcf.tsv
paste -d"\t" $INTERSECT_OUT/deboever_header_bed.tsv $INTERSECT_OUT/deboever_header_vcf.tsv > $INTERSECT_OUT/deboever_combined_header.tsv

cat $INTERSECT_OUT/deboever_combined_header.tsv > $INTERSECT_OUT/deboever_imputed_overlapping_filtered_header.bed
cat $INTERSECT_OUT/deboever_imputed_overlapping_filtered.bed >> $INTERSECT_OUT/deboever_imputed_overlapping_filtered_header.bed


### Use this file to test for trends in gene
awk 'BEGIN{FS=OFS="\t"}{print($31, $15, $25)}' $INTERSECT_OUT/deboever_imputed_overlapping_filtered_header.bed | sed 's/\..*\t/\t/g' > $INTERSECT_OUT/deboever_gene_snp_list.tsv



