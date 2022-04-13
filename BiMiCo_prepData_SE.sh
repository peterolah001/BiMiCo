#!/bin/bash

## Prepare Fastq files for BiMiCo pipeline - QC and filter

## ARG1 - exact file extension of READ1 files to process in drectory, e.g. .fastq (no quotes)

ARG1=$1

 for infile in *$1; do

   infile2=${infile/$1/$2}
   outfile1=${infile/$1/_filt$1}
   htmlrep=${infile/$1/.fastp.html}
   jsonrep=${infile/$1/.fastp.json}

fastp \
   -w 4 \
   -i $infile \
   -o $outfile1 \
   --dont_overwrite \
   -h $htmlrep \
   -j $jsonrep ;

done
