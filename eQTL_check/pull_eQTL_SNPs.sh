#!/bin/bash


vcf=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/SNPgenotypes/merged_imputed_AllChrs_iPSC.vcf
kilpinen=/directflow/SCCGGroupShare/projects/DrewNeavin/References/iPSC_eQTLs/KilpineniPSCeQTLs.txt
deboever=/directflow/SCCGGroupShare/projects/DrewNeavin/References/iPSC_eQTLs/DeBoeveriPSCeQTLs.txt

OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/References/iPSC_eQTLs
OUT_TMP=/directflow/SCCGGroupShare/projects/DrewNeavin/References/iPSC_eQTLs/tmp

mkdir -p $OUT_TMP

INTERSECT_OUT=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check/KilpinenOverlap 

mkdir -p $INTERSECT_OUT


##### KILPINEN DATA #####
#### Make bed file ####

### Get just first column ###
awk 'BEGIN{OFS="\t"}{print($1, $2 -1, $2)}' $kilpinen | \
    sed 's/chr\t-1\tpos/#chrom\tstart\tend/g' > $OUT_TMP/kilpinen_tmp.tsv



### Add back rest of the columns ###
cut -f 3- $kilpinen > $OUT_TMP/other_columns.tsv


paste $OUT_TMP/kilpinen_tmp.tsv $OUT_TMP/other_columns.tsv > $OUT/kilpinen.bed


### Remove unacceptable records (copy njmber of inversions) from vcf ###


sed '/<CN/d' $vcf > $INTERSECT_OUT/merged_imputed_AllChrs_iPSC_.vcf
sed -i '/<INV>/d' $INTERSECT_OUT/merged_imputed_AllChrs_iPSC_.vcf


### Intersect with unfiltered vcf ###
conda activate bedtools

    bedtools intersect -a $OUT/kilpinen.bed -b $INTERSECT_OUT/merged_imputed_AllChrs_iPSC_.vcf -wa -wb > $INTERSECT_OUT/kilpinen_imputed_overlapping.bed

conda deactivate


### Filter on individuals and MAF for our lines ###
VCF_R2=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/SNPgenotypes/merged_imputed_AllChrs_iPSC_R2_0.3.vcf

sed '/<CN/d' $VCF_R2 > $INTERSECT_OUT/merged_imputed_AllChrs_iPSC_R2_0.3.vcf
sed -i '/<INV>/d' $INTERSECT_OUT/merged_imputed_AllChrs_iPSC_R2_0.3.vcf


INDIVS=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/SNPgenotypes/lines.txt

conda activate vcftools

    vcftools --vcf $INTERSECT_OUT/merged_imputed_AllChrs_iPSC_R2_0.3.vcf --keep $INDIVS --recode --recode-INFO-all --mac 1 -o $INTERSECT_OUT/merged_imputed_AllChrs_iPSC_R2_0.3_filtered

conda deactivate

bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES || COUNT(GT="RA")=N_SAMPLES' $INTERSECT_OUT/merged_imputed_AllChrs_iPSC_R2_0.3_filtered.recode.vcf > $INTERSECT_OUT/merged_imputed_AllChrs_iPSC_R2_0.3_filtered_diff_genotypes.vcf


### Redo intersection with these SNPs ##$#
conda activate bedtools

    bedtools intersect -a $OUT/kilpinen.bed -b $INTERSECT_OUT/merged_imputed_AllChrs_iPSC_R2_0.3_filtered_diff_genotypes.vcf -wa -wb > $INTERSECT_OUT/kilpinen_imputed_overlapping_filtered.bed

conda deactivate

### Use this file to test for trends in gene
awk 'BEGIN{FS=OFS="\t"}{print($24, $7)}' $INTERSECT_OUT/kilpinen_imputed_overlapping_filtered.bed | sed 's/\..*//g' > $INTERSECT_OUT/gene_snp_list.tsv

sed -i '1i snp\tgene' $INTERSECT_OUT/gene_snp_list.tsv







##### DEBOEVER DATA #####
#### Make bed file ####

### Get just first three column for bed ###
awk 'BEGIN{OFS="\t"}{print($2, $3, $4)}' $deboever | \
    sed 's/^chr//g' |
    sed 's/start\tend\tmarker_id/chrom\tstart\tend/g' > $OUT_TMP/deboever_tmp.tsv



### Add back rest of the columns ###
cut -f 5- $deboever > $OUT_TMP/deboever_other_columns.tsv


paste $OUT_TMP/deboever_tmp.tsv $OUT_TMP/deboever_other_columns.tsv > $OUT/deboever.bed


### Intersect with these SNPs ##$#
conda activate bedtools

    bedtools intersect -a $OUT/deboever.bed -b $INTERSECT_OUT/merged_imputed_AllChrs_iPSC_R2_0.3_filtered_diff_genotypes.vcf -wa -wb > $INTERSECT_OUT/deboever_imputed_overlapping_filtered.bed

conda deactivate


head -n 1 $OUT/deboever.bed > $INTERSECT_OUT/deboever_header_bed.tsv
grep "#CHROM" $INTERSECT_OUT/merged_imputed_AllChrs_iPSC_R2_0.3_filtered_diff_genotypes.vcf > $INTERSECT_OUT/deboever_header_vcf.tsv
paste -d"\t" $INTERSECT_OUT/deboever_header_bed.tsv $INTERSECT_OUT/deboever_header_vcf.tsv > $INTERSECT_OUT/deboever_combined_header.tsv

cat $INTERSECT_OUT/deboever_combined_header.tsv > $INTERSECT_OUT/deboever_imputed_overlapping_filtered_header.bed
cat $INTERSECT_OUT/deboever_imputed_overlapping_filtered.bed >> $INTERSECT_OUT/deboever_imputed_overlapping_filtered_header.bed


### Use this file to test for trends in gene
awk 'BEGIN{FS=OFS="\t"}{print($31, $15, $25)}' $INTERSECT_OUT/deboever_imputed_overlapping_filtered_header.bed | sed 's/\..*\t/\t/g' > $INTERSECT_OUT/deboever_gene_snp_list.tsv



