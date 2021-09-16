#!/bin/bash


# if [[ $pool == "DRENEA_1" ]] || [[ $pool == "DRENEA_2" ]] || [[ $pool == "DRENEA_3" ]] || [[ $pool == "DRENEA_4" ]] || [[ $pool == "DRENEA_5" ]] || [[ $pool == "DRENEA_6" ]] 
# then
# 	bamtofastq-1.3.2 --reads-per-fastq=500000000 $TENxDIR/$pool/outs/possorted_genome_bam.bam $OUT/$pool/possorted_genome_bam
# else
# 	bamtofastq-1.3.2 --reads-per-fastq=500000000 $TENxDIR/$pool/$pool/outs/possorted_genome_bam.bam $OUT/$pool/possorted_genome_bam
# fi


cat $OUT/$pool/possorted_genome_bam/*/*I1* > $OUT/$pool/$pool\_L001_I1_001.fastq.gz
cat $OUT/$pool/possorted_genome_bam/*/*R1* > $OUT/$pool/$pool\_L001_R1_001.fastq.gz
cat $OUT/$pool/possorted_genome_bam/*/*R2* > $OUT/$pool/$pool\_L001_R2_001.fastq.gz