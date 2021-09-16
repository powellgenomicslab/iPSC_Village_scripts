


SIDX="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/preprocess/gencode.v38.annotation.expanded.sidx"
tx2gene="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/preprocess/gencode.v38.annotation.expanded.tx2gene.tsv"


eval "$(conda shell.bash hook)"
conda activate salmon

	salmon alevin -l ISR -i $SIDX \
		-1 $OUT/$pool/$pool\_L001_R1_001.fastq.gz \
		-2 $OUT/$pool/$pool\_L001_R2_001.fastq.gz \
		-o $OUT/$pool/alevin_out --tgMap $tx2gene \
		--chromium --dumpFeatures --expectCells 20000 -p $T

conda deactivate