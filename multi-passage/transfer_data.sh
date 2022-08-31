

INDIR=/directflow/SCCGGroupShare/projects/himaro/demuxafy/demultiplex_scRNA-seq/log_new/combined_results/AtLeastHalfSinglet
OUTDIR=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/demultiplexed

mkdir -p $OUTDIR

# for VIL in `ls $INDIR`
# do
#     cp -R $VIL $OUTDIR
# done

### Asked Himanshi to redo the overlap with demuxlet included since it allowed us to identify additionao 
mkdir -p $OUTDIR/with_demuxlet

for VIL in `ls $INDIR`
do
    cp -R $INDIR/$VIL $OUTDIR/with_demuxlet
done



conda activate baseR402

for village in `ls /directflow/SCCGGroupShare/projects/himaro/demuxafy/demultiplex_scRNA-seq/log_new/demuxlet_updated_barcodes/`
do
    mkdir -p /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/demultiplexed/with_demuxlet/updated_2022_06_26/atleasthalf_singlet/$village
    mkdir -p /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/demultiplexed/with_demuxlet/updated_2022_06_26/majority_singlet/$village

    Rscript /directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Demultiplexing_Doublet_Detecting_Docs/scripts/Combine_Results.R \
        -o /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/demultiplexed/with_demuxlet/updated_2022_06_26/majority_singlet/$village/majority_singlet_new.tsv \
        -d /directflow/SCCGGroupShare/projects/himaro/demuxafy/demultiplex_scRNA-seq/log_new/demuxlet_updated_barcodes/$village \
        -f /directflow/SCCGGroupShare/projects/himaro/demuxafy/demultiplex_scRNA-seq/log_new/freemuxlet_new_barcode/$village \
        -u /directflow/SCCGGroupShare/projects/himaro/demuxafy/demultiplex_scRNA-seq/log_new/souporcell_new_barcode/$village \
        -v /directflow/SCCGGroupShare/projects/himaro/demuxafy/demultiplex_scRNA-seq/log_new/vireo_new_barcode/$village \
        --method "MajoritySinglet"

        Rscript /directflow/SCCGGroupShare/projects/DrewNeavin/Demultiplex_Benchmark/Demultiplexing_Doublet_Detecting_Docs/scripts/Combine_Results.R \
        -o /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/demultiplexed/with_demuxlet/updated_2022_06_26/atleasthalf_singlet/$village/atleasthalf_singlet_new.tsv \
        -d /directflow/SCCGGroupShare/projects/himaro/demuxafy/demultiplex_scRNA-seq/log_new/demuxlet_updated_barcodes/$village \
        -f /directflow/SCCGGroupShare/projects/himaro/demuxafy/demultiplex_scRNA-seq/log_new/freemuxlet_new_barcode/$village \
        -u /directflow/SCCGGroupShare/projects/himaro/demuxafy/demultiplex_scRNA-seq/log_new/souporcell_new_barcode/$village \
        -v /directflow/SCCGGroupShare/projects/himaro/demuxafy/demultiplex_scRNA-seq/log_new/vireo_new_barcode/$village \
        --method "AtLeastHalfSinglet"
done