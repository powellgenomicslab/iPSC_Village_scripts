

DIR=/directflow/SCCGGroupShare/projects/himaro/demuxafy/demultiplex_scRNA-seq/log_dir/vireo_new_barcode
SIF=/directflow/SCCGGroupShare/projects/himaro/demuxafy/image/Demuxafy.sif
PIPELINE="/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/scripts/multi-passage/test_vireo_ambient/vireo.sh"

THREADS=4
N=18


for pool in `ls /directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/multi-passage/resequencing/220405/`
do

    FIELD="GP"


    VIREO_INDIR=/directflow/SCCGGroupShare/projects/himaro/demuxafy/demultiplex_scRNA-seq/log_new/vireo_new_barcode/$pool
    VIREO_OUTDIR=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/multi-passage/vireo_test/$pool
    LOG=$VIREO_OUTDIR/logs
    mkdir -p $LOG



    qsub -S /bin/bash \
        -cwd \
        -pe smp 4 \
        -N vireo \
        -q short.q \
        -l mem_requested=8G \
        -l tmp_requested=8G \
        -e $LOG \
        -o $LOG \
        -r yes \
        -j y \
        -M d.neavin@garvan.org.au \
        -v DIR=$DIR,SIF=$SIF,VIREO_INDIR=$VIREO_INDIR,VIREO_OUTDIR=$VIREO_OUTDIR,FIELD=$FIELD \
        -V \
        -C '' $PIPELINE


done