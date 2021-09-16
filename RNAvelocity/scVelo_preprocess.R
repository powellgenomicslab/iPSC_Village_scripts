## Date: 28 May, 2021
## Author: Drew Neavin
## Reason: RNA velocity on iPSC village samples for pseudotime following the Allevin tutorial https://combine-lab.github.io/alevin-tutorial/2020/alevin-velocity/

library(Biostrings)
library(BSgenome)
library(eisaR)
library(GenomicFeatures)
library(SummarizedExperiment)
library(tximeta)
library(rjson)
library(reticulate)
library(SingleCellExperiment)
library(scater)


outdir <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/output/scVelo/preprocess/"
dir.create(outdir, recursive = TRUE)

gtf <- "/directflow/SCCGGroupShare/projects/DrewNeavin/References/GenCode/GRCh38/gencode.v38.annotation.gtf.gz"
fasta <- "/directflow/SCCGGroupShare/projects/DrewNeavin/References/GenCode/GRCh38/GRCh38.primary_assembly.genome.fa"
pool_dir_10x <- "/directflow/SCCGGroupShare/projects/DrewNeavin/iPSC_Village/data/Expression_200128_A00152_0196_BH3HNFDSXY/GE/"


##################################################
##### Step 1. Generate reference fasta files #####
##################################################
##### Read in gtf file as GRanges Object #####
grl <- eisaR::getFeatureRanges(
  gtf = gtf,
  featureType = c("spliced", "intron"), 
  intronType = "separate", 
  flankLength = 90L, 
  joinOverlappingIntrons = FALSE, 
  verbose = TRUE
)


##### Make fasta file overlapping gtf #####
genome <- Biostrings::readDNAStringSet(
    fasta
)

names(genome) <- sapply(strsplit(names(genome), " "), .subset, 1)


seqs <- GenomicFeatures::extractTranscriptSeqs(
  x = genome, 
  transcripts = grl
)


Biostrings::writeXStringSet(
    seqs, filepath = paste0(outdir, "gencode.v38.annotation.expanded.fa")
)


##### Write Explanded annotation file #####
eisaR::exportToGtf(
  grl, 
  filepath = "gencode.v38.annotation.expanded.gtf"
)


##### Write spliced and unspliced seqwuences #####
head(metadata(grl)$corrgene)


write.table(
    metadata(grl)$corrgene, 
    file = paste0(outdir,"gencode.v38.annotation.expanded.features.tsv"),
    row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t"
)


##### Write file for mapping transcripts to genes #####
df <- eisaR::getTx2Gene(
    grl, filepath = paste0(outdir,"gencode.v38.annotation.expanded.tx2gene.tsv")
)


################################################
##### Step 2. Index the reference features #####
################################################
system(paste0('grep ">" ',fasta,' | cut -d ">" -f 2 | cut -d " " -f 1 > ', outdir, 'GRCh38.primary_assembly.genome.chrnames.txt'))

system(paste0('cat ',outdir, 'gencode.v38.annotation.expanded.fa ', fasta, '> ', outdir, 'combined_fasta.fa'))

### Run salmon index using submission system (qsub): salmon_index.qsub


##### Make tximeta linked transcriptome #####
tximeta::makeLinkedTxome(
  indexDir = paste0(outdir,"gencode.v38.annotation.expanded.sidx"), 
  source = "GENCODE", genome = "GRCh38", 
  organism = "Homo sapiens", release = "38", 
  fasta = paste0(outdir,"gencode.v38.annotation.expanded.fa"), 
  gtf = paste0(outdir,"gencode.v38.annotation.expanded.gtf"), 
  write = TRUE, jsonFile = "gencode.v38.annotation.expanded.json")


rjson::fromJSON(file = "gencode.v38.annotation.expanded.json")

########################################
##### Step 3. Quantify with alevin #####
########################################
### Run with bam2fastq.sh for each pool
