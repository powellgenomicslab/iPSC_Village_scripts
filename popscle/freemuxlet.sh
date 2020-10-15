#!/bin/bash

### If running on wolfpack7
eval "$(conda shell.bash hook)"
conda activate popscle

SAMPLE_INFO="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/SampleInfo.txt"
NAMES=$(cat $SAMPLE_INFO|cut -f1|tr '\n' ':')
NUMBERS=$(cat $SAMPLE_INFO|cut -f3|tr '\n' ':')
SAMPLE=$(echo $NAMES | cut -d: -f$SGE_TASK_ID)
N=$(echo $NUMBERS | cut -d: -f$SGE_TASK_ID)

BARCODES="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/200128_A00152_0196_BH3HNFDSXY/GE/$SAMPLE/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"

mkdir -p $OUT/$SAMPLE

echo $SAMPLE
echo $N

if [[ $stage == "pileup" ]]
then
    VCF=/directflow/SCCGGroupShare/projects/DrewNeavin/References/GRCh38SNPvcfs1000genomes/MergedAutosomesFiltered.recode.SortedReorderedMAF0.1.vcf # Exons only
    popscle dsc-pileup --sam $DATA/$SAMPLE/outs/possorted_genome_bam.bam --vcf $VCF --group-list $BARCODES --out $OUT/$SAMPLE/dsc.pileup0.1

elif [[ $stage == "freemuxlet" ]]
then
    echo "Freemuxlet with MAF 0.1"
    popscle freemuxlet --plp $DIR/output/popscle/pileup/$SAMPLE/dsc.pileup0.1 --out $OUT/$SAMPLE/MAF0.1freemuxletOUT --group-list $BARCODES --nsample $N

elif [[ $stage == "demuxlet" ]]
then
    # ### For with imputation at MAF 0.1, GP field
    # /share/ClusterShare/software/contrib/sccg/anaconda3/envs/demuxlet/bin/demuxlet --vcf $DIR/data/SNPgenotypes/merged_imputed_AllChrs_iPSC_R2_0.3_exons_filteredLines_hg38_MAF0.1.vcf --sam $DATA/$SAMPLE/outs/possorted_genome_bam.bam --field GP --group-list $BARCODES --out $OUT/$SAMPLE/ImputedMAF0.1_GPdemuxletOUT

    ### For with imputation no MAF, GP field
    /share/ClusterShare/software/contrib/sccg/anaconda3/envs/demuxlet/bin/demuxlet --vcf $DIR/data/SNPgenotypes/merged_imputed_AllChrs_iPSC_R2_0.3_exons_filteredLines_hg38.vcf --sam $DATA/$SAMPLE/outs/possorted_genome_bam.bam --field GP --group-list $BARCODES --out $OUT/$SAMPLE/Imputed_GPdemuxletOUT


    ### Note: this doesn't take very much memory, less than 1G of mem
fi