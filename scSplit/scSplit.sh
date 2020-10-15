#!/bin/bash

echo $T


export PYTHONPATH="/directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/scSplit_new/lib/python3.7/site-packages:$PYTHONPATH"

FASTA="/directflow/SCCGGroupShare/projects/DrewNeavin/References/ENSEMBLfasta/GRCh38/genome.fa"
FAI="/directflow/SCCGGroupShare/projects/DrewNeavin/References/ENSEMBLfasta/GRCh38/genome.fa.fai"
# scSplit="/home/drenea/.conda/envs/scSplit/lib/python3.7/site-packages/scSplit/scSplit"
# scSplit="/share/ClusterShare/software/contrib/sccg/miniconda3/envs/scSplit_w7/lib/python3.7/site-packages/scSplit/scSplit"
# scSplit="/share/ClusterShare/software/contrib/sccg/miniconda3/envs/scSplit_new/lib/python3.7/site-packages/scSplit/scSplit"
scSplit="/share/ScratchGeneral/drenea/tools/scSplit/scSplit"

N=`head -n $SGE_TASK_ID /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/SampleInfo.txt | tail -n 1 | awk 'BEGIN{FS="\t"}{print($3)}'`
echo $N
out=`head -n $SGE_TASK_ID /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/SampleInfo.txt | tail -n 1 | awk 'BEGIN{FS="\t"}{print($1)}'`
echo $out
BAM="$DATA/$out/outs/possorted_genome_bam.bam"
echo $BAM
BARCODES="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/200128_A00152_0196_BH3HNFDSXY/GE/$out/outs/filtered_feature_bc_matrix/barcodes.tsv"
echo $BARCODES
SNVs="/directflow/SCCGGroupShare/projects/DrewNeavin/References/scSplit_common_SNVs/common_snvs_hg38"


if [[ $stage == "prepare_files" ]]
then

    eval "$(conda shell.bash hook)"
    conda activate scSplit_new

    echo "Preparing Files"
    mkdir -p $OUT/InputFiles/$out
    mkdir -p $OUT/$out

    # ## Convert bam into sam
    # samtools view -h $BAM > $OUT/InputFiles/$out/SAM_body.sam
    # echo "The size pre-filtering"
    # wc -l $OUT/InputFiles/$out/SAM_body.sam

    # ## Remove sam files (HUGE!!!)
    # rm $OUT/InputFiles/$out/SAM_body.sam

    ## Save the header lines
    samtools view -@ $T -H $BAM > $OUT/InputFiles/$out/SAM_header

    ## Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
    samtools view -@ $T -S -q 10 -F 3844 $BAM | LC_ALL=C grep -F -f $BARCODES > $OUT/InputFiles/$out/filtered_SAM_body

    ## Combine header and body
    cat $OUT/InputFiles/$out/SAM_header $OUT/InputFiles/$out/filtered_SAM_body > $OUT/InputFiles/$out/filtered.sam

    ## Convert filtered.sam to BAM format
    samtools view -@ $T -b $OUT/InputFiles/$out/filtered.sam > $OUT/InputFiles/$out/filtered.bam

    # ## Remove sam files (HUGE!!!)
    rm $OUT/InputFiles/$out/filtered_SAM_body
    rm $OUT/InputFiles/$out/SAM_header
    rm $OUT/InputFiles/$out/filtered.sam

    ## Remove duplicates
    samtools rmdup $OUT/InputFiles/$out/filtered.bam $OUT/InputFiles/$out/dedup_filtered.bam

    ## Sort
    samtools sort -@ $T -o $OUT/InputFiles/$out/possort_dedup_filtered.bam $OUT/InputFiles/$out/dedup_filtered.bam

    ## Generate index for BAM
    samtools index $OUT/InputFiles/$out/possort_dedup_filtered.bam

    ## Remove the extra bam files
    rm $OUT/InputFiles/$out/dedup_filtered.bam
    rm $OUT/InputFiles/$out/filtered.bam

elif [[ $stage == "freebayes" ]]
then

    eval "$(conda shell.bash hook)"
    conda activate freebayes

    # echo "Starting Freebayes"
    # ## Use freebayes to identify SNV
    # fasta_generate_regions.py $FAI 100000 > $OUT/InputFiles/$out/regions_file
    # freebayes-parallel $OUT/InputFiles/$out/regions_file 40 -f $FASTA -iXu -C 20 -q 3 $OUT/InputFiles/$out/possort_dedup_filtered.bam > $OUT/InputFiles/$out/var_count20.vcf
    

    # ## Get just the high quality calls
    # vcftools --gzvcf $OUT/InputFiles/$out/var_count20.vcf --minQ 30 --recode --recode-INFO-all --out $OUT/InputFiles/$out/var_count20_qual30.vcf


    ## Pull just the SNPs with MAF between 0.1 and 0.9
    bcftools filter --include 'AF<=0.9 && AF>=0.1' -O v --output $OUT/InputFiles/$out/var_count20_qual30_filtered.vcf $OUT/InputFiles/$out/var_count20_qual30.vcf.recode.vcf

elif [[ $stage == "AlleleMatrices" ]]
then

    eval "$(conda shell.bash hook)"
    conda activate scSplit_new

    echo "Creating Allele Matrices"
    ## Run python script "matrices.py" and get two .csv files ("ref_filtered.csv" and "alt_filtered.csv") as output.
    cd $OUT/$out

    ### When used freebayes with counts 20
    python $scSplit count -c $SNVs -v $OUT/InputFiles/$out/var_count20_qual30_filtered.vcf -i $OUT/InputFiles/$out/possort_dedup_filtered.bam -b $BARCODES -r $OUT/$out/ref_filtered_count20.csv -a $OUT/$out/alt_filtered_count20.csv 


elif [[ $stage == "Demultiplex" ]]
then
    echo "Starting demultiplexing"
    ## Use the two generated allele counts matrices files to demultiplex the cells into different samples. Doublet sample will not have the same sample ID every time, which will be explicitly indicated in the log file
    mkdir -p $OUT/$out/Count20
    cd $OUT/$out/Count20
    python $scSplit run -r $OUT/$out/ref_filtered_count20.csv -a $OUT/$out/alt_filtered_count20.csv -n $N
    
    # mkdir -p $OUT/$out/Map2Genotypes/
    # cd $OUT/$out/Map2Genotypes/
    # python $scSplit run -v $DIR/output/vireo/Imputed/$SAMPLE/SUBMergedMAF0.05.dose.vcf.gz -r $OUT/$out/ref_filtered.csv -a $OUT/$out/alt_filtered.csv -n $N
    ### Does not work -- assigns as donor 1, 2 etc.


elif [[ $stage == "scSplitGenotypes" ]]
then
    echo "Getting Genotypes"
    ### Generate genotypes based on the split result ###
    ## Count 20 ##
    cd $OUT/$out/Count20
    python $scSplit genotype -r $OUT/$out/ref_filtered_count20.csv -a $OUT/$out/alt_filtered_count20.csv -p $OUT/$out/Count20/scSplit_P_s_c.csv

fi