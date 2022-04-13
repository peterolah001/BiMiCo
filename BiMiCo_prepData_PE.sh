#!/bin/bash

## Prepare Fastq files for BiMiCo pipeline - QC and filter

## ARG1 - exact file extension of READ1 files to process in drectory, e.g. .fastq (no quotes)
## ARG2 - exact file extension of READ2 files to process in drectory, e.g. .fastq (no quotes)

ARG1=$1
ARG2=$2

 for infile in *$1; do

   infile2=${infile/$1/$2}
   outfile1=${infile/$1/_filt$1}
   outfile2=${infile/$1/_filt$2}
   htmlrep=${infile/$1/.report.html}
   jsonrep=${infile/$1/.report.json}

fastp \
   -w 4 \
   -i $infile \
   -I $infile2 \
   -o ./Trimmed_reads/$outfile1 \
   -O ./Trimmed_reads/$outfile2 \
   --dont_overwrite \
   -h $htmlrep \
   -j $jsonrep
