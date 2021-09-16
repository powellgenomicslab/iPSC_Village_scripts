Kilpinen="/directflow/SCCGGroupShare/projects/DrewNeavin/References/iPSC_eQTLs/KilpineniPSCeQTLs.txt"
DeBoever="/directflow/SCCGGroupShare/projects/DrewNeavin/References/iPSC_eQTLs/DeBoeveriPSCeQTLs.txt"
OUT="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check"
VCF="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/SNPgenotypes/merged_imputed_AllChrs_iPSC_R2_0.3.vcf"

mkdir -p $OUT

gene="CHCHD2"
ensg="ENSG00000106153"

head -n 1 $DeBoever > $OUT/DeBoever_CHCHD2.tsv
grep $gene $DeBoever >> $OUT/DeBoever_CHCHD2.tsv
awk 'BEGIN{OFS=""}{print($2,":",$3,"-",$4)}' $OUT/DeBoever_CHCHD2.tsv > $OUT/DeBoever_CHCHD2_4liftover.tsv


head -n 1 $Kilpinen > $OUT/Kilpinen_CHCHD2.tsv
grep $ensg $Kilpinen >> $OUT/Kilpinen_CHCHD2.tsv
awk 'BEGIN{FS=OFS="\t"}{print($1,$2-1,$2)}' $OUT/Kilpinen_CHCHD2.tsv > $OUT/Kilpinen_CHCHD2_4liftover.tsv


### Lifted over using UCSC lifteover browser


### Convert to no chr and second number to grep from vcf
##  DeBoever
sed -i 's/:/\t/g' $OUT/DeBoever_CHCHD2_lifted_hg38.bed
sed -i 's/-/\t/g' $OUT/DeBoever_CHCHD2_lifted_hg38.bed
awk 'BEGIN{FS=OFS="\t"}{print($1, $3)}' $OUT/DeBoever_CHCHD2_lifted_hg38.bed > $OUT/DeBoever_CHCHD2_lifted_hg38.tsv
sed -i 's/chr//g' $OUT/DeBoever_CHCHD2_lifted_hg38.tsv

## Kilpinen
awk 'BEGIN{FS=OFS="\t"}{print($1, $3)}' $OUT/Kilpinen_CHCHD2_lifted_hg38.bed > $OUT/Kilpinen_CHCHD2_lifted_hg38.tsv
sed -i 's/chr//g' $OUT/Kilpinen_CHCHD2_lifted_hg38.tsv

grep "#CHROM" $VCF > $OUT/eQTL_SNP_genotypes.vcf
grep -f $OUT/DeBoever_CHCHD2_lifted_hg38.tsv $VCF >> $OUT/eQTL_SNP_genotypes.vcf
grep "56018361" $VCF

### The Kilpinen SNP does not exist in our vcf - find SNPs in LD and pull those from VCF
## Use snpsnap to find snps in LD, then filter for chromosome 7
## 1. Pull the new SNP from column 3 with awk
awk 'BEGIN{FS=OFS="\t"}{print($3)}' $OUT/snpsnap_r2_results.tsv | sed 's/:/\t/g' | awk 'BEGIN{FS=OFS="\t"}{print($1, $2-1, $2)}' | grep "^7	" | sort -u > $OUT/snpsnap_locations.bed


## 2. Make Kilpinen into bed file
awk 'BEGIN{FS=OFS="\t"}{print($1, $2-1, $2)}' $OUT/Kilpinen_CHCHD2_lifted_hg38.tsv > $OUT/Kilpinen_CHCHD2_lifted_hg38.bed

## 3. use bedtools to pull just the snps within 1Mb of main SNP
conda activate bedtools
bedtools window -a $OUT/Kilpinen_CHCHD2_lifted_hg38.bed -b $OUT/snpsnap_locations.bed -w 500000 | awk 'BEGIN{FS=OFS="\t"}{print($4,$5,$6)}' | sed 's/^/chr/g' > $OUT/SNPs_in_LD_Kilpinen.bed
conda deactivate

## Do liftover to hg38 from hg19 on ucsc 

## 4. Get just the chr and bp from the file (columns 1 and 3)
awk 'BEGIN{FS=OFS="\t"}{print($1,$3)}' $OUT/SNPs_in_LD_Kilpinen_hg38.bed | sed 's/chr//g' > $OUT/SNPs_in_LD_Kilpinen.tsv


## Grep for those SNPs in vcf
grep -f $OUT/SNPs_in_LD_Kilpinen.tsv $VCF > $OUT/eQTL_SNP_genotypes_Kilpinen.vcf
grep "^7	5597" $VCF | head
