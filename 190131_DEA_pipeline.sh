#!/bin/bash

# Standard RNAseq data processing pipeline as of 190131
# fastqc --> trim_galore --> multiQC --> STAR alignment --> featureCounts

# parameters:
# param 1: root directory
# should contain raw-data subdirectory with gzipped fastq files
# param 2: stranding
# 0 = unstranded or paired, 1 = first-strand, 2 = second-strand
# param 3: number of processors
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

#first determine if the files are paired
cd raw_data
isPaired=0
pairedFiles=`find . -type f -name "*_2.fastq.gz"` #if there are paired files, they will have this format
if [ ${#pairedFiles} -gt 0 ]
then
	isPaired=1
	echo 'files are in paired-end format'
else
	echo 'files are in single-end format'
fi

# initial QC
cd "$wd"
mkdir -p fastqc
mkdir -p ./fastqc/before_trim
mkdir -p ./fastqc/after_trim
fastqc -t $THREADS -o ./fastqc/before_trim ./raw_data/*.fastq.gz
multiqc ./fastqc/before_trim/* -o ./fastqc/before_trim/

# trimming and final QC (using GNU parallel)
mkdir -p trimmed

if [ $isPaired -eq 0 ]
then
	parallel -j $THREADS trim_galore \
	--quality $QCUTOFF --phred33 -o ./trimmed \
	--fastqc -o ./fastqc/after_trim ::: ./raw_data/*.fastq.gz
	multiqc ./fastqc/after_trim/* -o ./fastqc/after_trim/
else
	#--xapply only runs pairs, as opposed to every combination of outputs
	parallel --xapply -j $THREADS trim_galore --paired \
	--quality $QCUTOFF --phred33 -o ./trimmed \
	--fastqc --fastqc_args "-outdir ./fastqc/after_trim" \
	::: ./raw_data/*_1.fastq.gz ::: ./raw_data/*_2.fastq.gz
	multiqc ./fastqc/after_trim/* -o ./fastqc/after_trim/
fi

# alignment
mkdir -p alignment
if [ $isPaired -eq 0 ]
then
	for file in ./raw_data/*trimmed.fq.gz
	do
	  prefix=$(basename -- "$file")
	  prefix="./alignment/""${prefix%_trimmed.fastq*}"
	  STAR --runMode alignReads --genomeDir $GENOME --runThreadN $THREADS \
	  --readFilesCommand 'gunzip -c' --readFilesIn $file \
	  --outFileNamePrefix $prefix
	done
else
	leftFiles=./raw_data/*_1_val_1.fq.gz
	for file in $leftFiles
	do
	  outPrefix=$(basename -- "$file")
	  inPrefix="${prefix%_1_val_1.fastq.gz}"
	  outPrefix="./alignment/""$outPrefix"
	  STAR --runMode alignReads --genomeDir $GENOME --runThreadN $THREADS \
	  --readFilesCommand 'gunzip -c' \
	  --readFilesIn "$inPrefix"_1_val_1.fastq.gz "$inPrefix"_2_val_2.fastq.gz \
	  --outFileNamePrefix $outPrefix
	done
fi

# read assigment / summarization
mkdir -p counts
Rscript ~/Documents/GitHub/morimotoLab/190131_Rsubread_pipeline.r \
"$WD" $STRANDING $THREADS "$GTF" $MINFRAG $COUNTMM "$OUTPUTPREFIX" \
$isPaired

# prepare for DE (allow user to do by hand for now)
mkdir -p DE_analysis
