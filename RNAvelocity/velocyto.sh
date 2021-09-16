#!/bin/bash



eval "$(conda shell.bash hook)"
conda activate /directflow/SCCGGroupShare/projects/walmus/.conda/envs/velocyto

if [[ $pool == "DRENEA_1" ]] || [[ $pool == "DRENEA_2" ]] || [[ $pool == "DRENEA_3" ]] || [[ $pool == "DRENEA_4" ]] || [[ $pool == "DRENEA_5" ]] || [[ $pool == "DRENEA_6" ]] 
then
	velocyto run10x --samtools-threads $T -vv $TENxDIR/$pool $GTF
else
	velocyto run10x --samtools-threads $T -vv $TENxDIR/$pool/$pool $GTF
fi

conda deactivate