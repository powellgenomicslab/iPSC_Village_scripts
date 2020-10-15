#!/bin/bash

echo $SGE_TASK_ID

SAMPLE_INFO="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/SampleInfo.txt"
SNPs="/directflow/SCCGGroupShare/projects/DrewNeavin/References/souporcell_commonSNPs/filtered_2p_1kgenomes_GRCh38.vcf"
FASTA="/directflow/SCCGGroupShare/projects/DrewNeavin/References/ENSEMBLfasta/GRCh38/genome.fa"
NAMES=$(cat $SAMPLE_INFO|cut -f1|tr '\n' ':')
pool=$(echo $NAMES | cut -d: -f$SGE_TASK_ID)
NUMBERS=$(cat $SAMPLE_INFO|cut -f3|tr '\n' ':')
N=$(echo $NUMBERS | cut -d: -f$SGE_TASK_ID)
BAM="$DATA/$pool/outs/possorted_genome_bam.bam"
OUT="$OUT/$pool"
BARCODES="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/200128_A00152_0196_BH3HNFDSXY/GE/$pool/outs/filtered_feature_bc_matrix/barcodes.tsv"
echo $BARCODES


mkdir -p $OUT

cd /directflow/SCCGGroupShare/projects/DrewNeavin/

if [[ $stage == "run_souporcell" ]]
then
    singularity exec /directflow/SCCGGroupShare/projects/DrewNeavin/souporcell.sif souporcell_pipeline.py -i $BAM -b $BARCODES -f $FASTA -t $T -o $OUT -k $N

elif [[ $stage == "fix_output_files" ]]
then
    sed 's/\t\t/\t/g' $OUT/clusters.tsv > $OUT/clusters_fixed_tabs.tsv
    sed -i 's/\t$//g' $OUT/clusters_fixed_tabs.tsv
fi