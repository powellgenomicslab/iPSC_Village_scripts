#!/bin/bash


VCF=/directflow/SCCGGroupShare/projects/DrewNeavin/References/GRCh38SNPvcfs1000genomes/MergedAutosomesFiltered.recode.SortedReorderedMAF0.1.vcf # Exons only

SAMPLE_INFO="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/SampleInfo.txt"
NAMES=$(cat $SAMPLE_INFO|cut -f1|tr '\n' ':')
out=$(echo $NAMES | cut -d: -f$SGE_TASK_ID)
NUMBERS=$(cat $SAMPLE_INFO|cut -f3|tr '\n' ':')
N=$(echo $NUMBERS | cut -d: -f$SGE_TASK_ID)
BAM="$DATA/$out/outs/possorted_genome_bam.bam"
echo $BAM
BARCODES="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/200128_A00152_0196_BH3HNFDSXY/GE/$out/outs/filtered_feature_bc_matrix/barcodes.tsv"
echo $BARCODES

if [[ $stage == "prepBAMS" ]]
then
    echo "Starting pre bam prep"
    # mkdir -p $OUT/$out/
    # cp $DATA/$out/outs/possorted_genome_bam.bam $OUT/$out/possorted_genome_bam_poolBarcode.bam
    # sed -i "s/\(CB:Z:.\{16\}-1\)/\1.$out/g" $OUT/$out/possorted_genome_bam_poolBarcode.bam
    # echo $OUT/$out/possorted_genome_bam_poolBarcode.bam >> $OUT/bamListFile.txt

elif [[ $stage == "pileup" ]]
then
    
    eval "$(conda shell.bash hook)"
    conda activate cellSNP

    echo "starting pileup"

    mkdir -p $OUT/$out
    cd $OUT/$out 

    cellSNP -s $BAM -b $BARCODES -o $OUT/$out/cellSNPpileupMAF0.1.vcf -R $VCF -p 20 --minMAF 0.1 --minCOUNT 20

elif [[ $stage == "vireoPrepareImputed" ]]
then
    eval "$(conda shell.bash hook)"
    conda activate vireo

    echo "starting new vireo preparation file generation with imputed SNPs"
    DONOR_GT_FILE="$DIR/data/SNPgenotypes/merged_imputed_AllChrs_iPSC_R2_0.3_exons_filteredLines_hg38_MAF0.1.vcf.gz"
    SAMPLE=$OUT/SampleIDlist/$out/sampleList.txt
    
    mkdir -p $OUT
    mkdir -p $OUT/$out

    bcftools view -R $DIR/output/vireo/$out/cellSNPpileupM0AF0.1.vcf.gz -Oz -o $OUT/$out/SUBMergedMAF0.1.dose.vcf.gz $DONOR_GT_FILE 


elif [[ $stage == "vireoGenotyped" ]]
then

    ps -o pid,ppid,lstart,etime,%cpu,%mem,rss,vsz,pmem,comm,maj_flt,min_flt

    eval "$(conda shell.bash hook)"
    conda activate vireo

    echo "starting new vireo"
    ##### When running with imputed SNPs 
    n_donor=`head -n $SGE_TASK_ID $SAMPLE_INFO | tail -n 1 | cut -f3`
    CELL_DATA="$OUT/$out/cellSNPpileupM0AF0.1.vcf.gz"
    DONOR_GT_FILE="$OUT/$out/SUBMergedMAF0.1.dose.vcf.gz"
    OUT=$OUT/$stage 
    OUT=$OUT/MAF0.1/$out 


    mkdir -p $OUT

    echo $out
    echo $n_donor
    echo $CELL_DATA
    echo $DONOR_GT_FILE

    vireo -c $CELL_DATA -d $DONOR_GT_FILE -o $OUT -t GP

    ps -o pid,ppid,lstart,etime,%cpu,%mem,rss,vsz,pmem,comm,maj_flt,min_flt

elif [[ $stage == "vireoUngenotyped" ]]
then
    echo "starting new vireo no genotype data"

    eval "$(conda shell.bash hook)"
    conda activate vireo

    ##### When running with MAF 0.1
    CELL_DATA="$OUT/$out/cellSNPpileupM0AF0.1.vcf.gz"
    n_donor=`head -n $SGE_TASK_ID $SAMPLE_INFO | tail -n 1 | cut -f3`
    OUT=$OUT/$stage/MAF0.1/$out
    mkdir -p $OUT

    echo $out
    echo $n_donor
    
    vireo -c $CELL_DATA -o $OUT -N $n_donor
fi

