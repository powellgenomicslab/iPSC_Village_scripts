#!/bin/bash

### Realized that called the wrong model and weren't really saving the correct model for qtl residuals
### Need to remove the files in icc and fit_models to force snakemake to rerun

DIR=/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/variance_partition_post_review/gene_separated

for file in `ls $DIR/residuals4qtl`
do
    base=`echo $file | sed 's/_residuals4qtl.rds//g'`
    echo $base
    rm $DIR/icc/$base\_icc.rds
    rm $DIR/fit_models/$base\_fitted_models.rds
done