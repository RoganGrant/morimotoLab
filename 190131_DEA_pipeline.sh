#!/bin/bash

# Standard RNAseq data processing pipeline as of 190131
# fastqc --> trim_galore --> multiQC --> STAR alignment --> featureCounts

# parameters:
# param 1: root directory
# should contain raw-data subdirectory with gzipped fastq files
# param 2: stranding
# 0 = unstranded or paired, 1 = first-strand, 2 = second-strand
# param 3 (optional): number of processors (defaults to 8)
# param 4: quality score cutoff for trimming
# param 5: directory of STAR-indexed genome to align to
# param 6: GTF file for read assignment
# param 7: minimum fragment length for featureCounts
# param 8: count multi-mapping reads (boolean integer)
# param 9: prefix to give output file (study name)

cd "$1" #set working directory
wd="$1"
STRANDNG=$2
THREADS=$3
QCUTOFF=$4
GENOME="$5"
GTF="$6"
MINFRAG=$7
COUNTMM=$8
OUTPUTPREFIX="$9"

# initial QC
mkdir -f fastqc
mkdir -f ./fastqc/before_trim
mkdir -f ./fastqc/after_trim
fastqc -t $THREADS -o ./fastqc/before_trim ./raw_data/*.fastq.gz
multiqc ./fastqc/before_trim/*

# trimming and final QC (using GNU parallel)
mkdir -f trimmed
parallel -j $THREADS trim_galore \
--quality $QSCORE --phred33 -o ./trimmed \
--fastqc -o ./fastqc/after_trim ::: ./raw_data/*.fastq.gz
multiqc ./fastqc/after_trim/*

# alignment
mkdir -f alignment
for file in ./raw_data/*trimmed.fq.gz
do
  prefix=$(basename -- "$i")
  prefix="./alingment/""${prefix%.fastq*}"
  STAR --runMode alignReads --genomeDir $GENOME --runThreadN $THREADS \
  --readFilesCommand 'gunzip -c' --readFilesIn $file \
  --outFileNamePrefix $prefix
done

# read assigment / summarization
mkdir -f counts
Rscript "~/Documents/GitHub/morimotoLab/190131_Rsubread_pipeline.r" \
"$WD" $STRANDING $THREADS "$GTF" $MINFRAG $COUNTMM "$OUTPUTPREFIX"

# prepare for DE (allow user to do by hand for now)
mkdir -f DE_analysis
