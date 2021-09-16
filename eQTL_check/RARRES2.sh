Kilpinen="/directflow/SCCGGroupShare/projects/DrewNeavin/References/iPSC_eQTLs/KilpineniPSCeQTLs.txt"
DeBoever="/directflow/SCCGGroupShare/projects/DrewNeavin/References/iPSC_eQTLs/DeBoeveriPSCeQTLs.txt"
OUT="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/eQTL_check"
VCF="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/SNPgenotypes/merged_imputed_AllChrs_iPSC_R2_0.3.vcf"

mkdir -p $OUT

gene="RARRES2"
ensg="ENSG00000106538"

PGx_SNP1="rs2159234"  #hg19 = 7	;
PGx_SNP2="rs3735175" #hg19 = 7	150020457;
PGx_SNP3="rs2108851" #hg19 = 7	150019634; hg38 = 7	150322545
PGx_SNP4="rs2108852" #hg19 = 7	150019728		;hg38 = 7	150322539


grep "#CHROM" $VCF | awk 'BEGIN{FS=OFS="\t"}{print($1,$2,$3,$4,$5,$6,$7,$8,$9,$23,$29,$21)}' > $OUT/$PGx_SNP1\_RARRES2.vcf
sed -n  "/7\t150020672/p" $VCF | awk 'BEGIN{FS=OFS="\t"}{print($1,$2,$3,$4,$5,$6,$7,$8,$9,$23,$29,$21)}' >> $OUT/$PGx_SNP1\_RARRES2.vcf

grep "#CHROM" $VCF | awk 'BEGIN{FS=OFS="\t"}{print($1,$2,$3,$4,$5,$6,$7,$8,$9,$23,$29,$21)}' > $OUT/$PGx_SNP2\_RARRES2.vcf
sed -n  "/7\t150020457/p" $VCF | awk 'BEGIN{FS=OFS="\t"}{print($1,$2,$3,$4,$5,$6,$7,$8,$9,$23,$29,$21)}' >> $OUT/$PGx_SNP2\_RARRES2.vcf #******** snp in ipscs avail FSA0006 = 0.001, MBE1006 = 1, TOB0421 = 1; C/T

grep "#CHROM" $VCF | awk 'BEGIN{FS=OFS="\t"}{print($1,$2,$3,$4,$5,$6,$7,$8,$9,$23,$29,$21)}' > $OUT/$PGx_SNP3\_RARRES2.vcf
sed -n  "/7\t150019634/p" $VCF | awk 'BEGIN{FS=OFS="\t"}{print($1,$2,$3,$4,$5,$6,$7,$8,$9,$23,$29,$21)}' >> $OUT/$PGx_SNP3\_RARRES2.vcf #******** snp in ipscs avail FSA0006 = 0.001, MBE1006 = 1, TOB0421 = 1; T/C

grep "#CHROM" $VCF | awk 'BEGIN{FS=OFS="\t"}{print($1,$2,$3,$4,$5,$6,$7,$8,$9,$23,$29,$21)}' > $OUT/$PGx_SNP4\_RARRES2.vcf
sed -n  "/7\t150019728/p" $VCF | awk 'BEGIN{FS=OFS="\t"}{print($1,$2,$3,$4,$5,$6,$7,$8,$9,$23,$29,$21)}' >> $OUT/$PGx_SNP4\_RARRES2.vcf


head -n 1 $DeBoever > $OUT/DeBoever_RARRES2.tsv
grep $gene $DeBoever >> $OUT/DeBoever_RARRES2.tsv
awk 'BEGIN{OFS="\t"}{print($2,$3)}' $OUT/DeBoever_RARRES2.tsv | sed 's/chr//g' > $OUT/DeBoever_RARRES2_4grep.tsv


head -n 1 $Kilpinen > $OUT/Kilpinen_RARRES2.tsv
grep $ensg $Kilpinen >> $OUT/Kilpinen_RARRES2.tsv
awk 'BEGIN{FS=OFS="\t"}{print($1,$2-1,$2)}' $OUT/Kilpinen_RARRES2.tsv > $OUT/Kilpinen_RARRES2_4liftover.tsv

grep "7	150002815" $VCF | awk 'BEGIN{FS=OFS="\t"}{print($1,$2,$3,$4,$5,$6,$7,$8,$9,$23,$29,$21)}' ## ***** snp in iPSCs FSA0006 = 0.055, MBE1006 = 0.943, TOB0421 = 0.982; T/A; beta = 0.454432062 in kilpinen


##### Get the 
grep "#CHROM" $VCF | awk 'BEGIN{FS=OFS="\t"}{print($1,$2,$3,$4,$5,$6,$7,$8,$9,$23,$29,$21)}' > $OUT/RARRES2_eQTL_SNP_genotypes.vcf
for snp in `cat $OUT/DeBoever_RARRES2_4grep.tsv`
do
	echo $snp
	# grep $snp $VCF | awk 'BEGIN{FS=OFS="\t"}{print($1,$2,$3,$4,$5,$6,$7,$8,$9,$23,$29,$21)}' >> $OUT/RARRES2_eQTL_SNP_genotypes.vcf
done



### Lifted over using UCSC lifteover browser from hg19 to hg38


### Convert to no chr and second number to grep from vcf
##  DeBoever
sed -i 's/:/\t/g' $OUT/DeBoever_RARRES2_lifted_hg38.bed
sed -i 's/-/\t/g' $OUT/DeBoever_RARRES2_lifted_hg38.bed
awk 'BEGIN{FS=OFS="\t"}{print($1, $3)}' $OUT/DeBoever_RARRES2_lifted_hg38.bed > $OUT/DeBoever_RARRES2_lifted_hg38.tsv
sed -i 's/chr//g' $OUT/DeBoever_RARRES2_lifted_hg38.tsv


grep "#CHROM" $VCF > $OUT/RARRES2_eQTL_SNP_genotypes.vcf
grep -f $OUT/DeBoever_RARRES2_lifted_hg38.tsv $VCF >> $OUT/RARRES2_eQTL_SNP_genotypes.vcf


## Kilpinen
awk 'BEGIN{FS=OFS="\t"}{print($1, $3)}' $OUT/Kilpinen_RARRES2_lifted_hg38.bed > $OUT/Kilpinen_RARRES2_lifted_hg38.tsv
sed -i 's/chr//g' $OUT/Kilpinen_RARRES2_lifted_hg38.tsv


### The Kilpinen SNP does not exist in our vcf - find SNPs in LD and pull those from VCF
## Use snpsnap to find snps in LD, then filter for chromosome 7
## 1. Pull the new SNP from column 3 with awk
awk 'BEGIN{FS=OFS="\t"}{print($3)}' $OUT/snpsnap_r2_results.tsv | sed 's/:/\t/g' | awk 'BEGIN{FS=OFS="\t"}{print($1, $2-1, $2)}' | grep "^7	" | sort -u > $OUT/snpsnap_locations.bed


## 2. Make Kilpinen into bed file
awk 'BEGIN{FS=OFS="\t"}{print($1, $2-1, $2)}' $OUT/Kilpinen_RARRES2_lifted_hg38.tsv > $OUT/Kilpinen_RARRES2_lifted_hg38.bed

## 3. use bedtools to pull just the snps within 1Mb of main SNP
conda activate bedtools
bedtools window -a $OUT/Kilpinen_RARRES2_lifted_hg38.bed -b $OUT/snpsnap_locations.bed -w 500000 | awk 'BEGIN{FS=OFS="\t"}{print($4,$5,$6)}' | sed 's/^/chr/g' > $OUT/SNPs_in_LD_Kilpinen.bed
conda deactivate

## Do liftover to hg38 from hg19 on ucsc 

## 4. Get just the chr and bp from the file (columns 1 and 3)
awk 'BEGIN{FS=OFS="\t"}{print($1,$3)}' $OUT/SNPs_in_LD_Kilpinen_hg38.bed | sed 's/chr//g' > $OUT/SNPs_in_LD_Kilpinen.tsv


## Grep for those SNPs in vcf
grep -f $OUT/SNPs_in_LD_Kilpinen.tsv $VCF > $OUT/eQTL_SNP_genotypes_Kilpinen.vcf
grep "^7	5597" $VCF | head
