---
<<<<<<< HEAD
# "Short tutorial for BiMiCo package basic functions"
=======
# Short tutorial for BiMiCo package basic functions
>>>>>>> 2aa0ee20b32363a8317aa9eb7e14d043abb24603
---



The BiMiCo (Biomap Microbiome Core-tools) package offers single-step preprocessing of 16S microbial marker gene sequencing raw data files for 454 and Illumina platforms. \
The current tutorial and example data is aimed at the simple preprocessing of 454-generated fastq data. For more options and details, please see the "BiMiCo_tutorial_LONG" document. (pls contact me for appropriate demo dataset)

# Install dependencies

- Install CRAN dependencies

`install.packages("BiocManager")` \
`install.packages("ggplot2")`

- Install Bioconductor dependencies

`BiocManager::install("Biobase")` \
`BiocManager::install("Biostrings")` \
`BiocManager::install("phyloseq")` \
`BiocManager::install("dada2")`

- Optional packages (required for "Long" tutorial/workflow)

`install.packages("fastqcr")` \
`install.packages("phangorn")` \
`BiocManager::install("decontam")`

# Load required packages

`library(BiMiCo)` \
`library(fastqcr)` \
`library(ggplot2)` \


# Specify INPUTS:

- fastq files (directory containing raw fastq files intended for analysis exclusively, no other files):

<<<<<<< HEAD
`rawfqs <- "./Demo_data_bimico/raw_fqs"`

- sample phenodata (sample names in rows, sample variables in columns):

`phedat <- read.csv("./Demo_data_bimico/phenodata_demo.csv", row.names=1)`
=======
`rawfqs <- "./Demo_data_bimico/raw_fqs"`

- sample phenodata (sample names in rows, sample variables in columns):

`phedat <- read.csv("./Demo_data_bimico/phenodata_demo.csv", row.names=1)`
>>>>>>> 2aa0ee20b32363a8317aa9eb7e14d043abb24603

- if not already contained in sample phenodata, please specify **PRIMER TYPE** (amplified region):

`phedat$primer_type <- "V3-V4"`

- taxonomy database to use with DADA2 (default is Greengenes v.13.5):

<<<<<<< HEAD
`txset <- "./Demo_data_bimico/gg_13_8_train_set_97.fa.gz"`
=======
`txset <- "./Demo_data_bimico/gg_13_8_train_set_97.fa.gz"`
>>>>>>> 2aa0ee20b32363a8317aa9eb7e14d043abb24603

# Specify OUTPUT folders:

- directory to write quality-filtered fastq files to:

<<<<<<< HEAD
`filtered_fqs <- "./Demo_data_bimico/filt_fqs"`

- directory to write various outputs(ASV table, figures & graphs) to:

`out_main <- "./Demo_data_bimico/results/"`
=======
`filtered_fqs <- "./Demo_data_bimico/filt_fqs"`

- directory to write various outputs(ASV table, figures & graphs) to:

`out_main <- "./Demo_data_bimico/results/"`
>>>>>>> 2aa0ee20b32363a8317aa9eb7e14d043abb24603


# Run single-command preprocessing for 454 reads. 
- *Tip #1: trimming the length of 454 sequencing reads is highly recommended due to the drop-off in quality; please set trim_read_length accordingly, see "Long" tutorial for more details on QC/trimming* \
- *Tip #2: Multithreading is enabled with the "mtthread=T" argument- not recommended when running from Rstudio.*

`toy_set <- bimico_454(rawfqs, txset, phedat, filtered_fqs, trim_read_length = 500, mtthread = F)`


# Generate a summary report based on the desired variable column of phenodata (case/control, sample source, sequencing batch...)

- UNDER DEVELOPMENT: specify phyloseq object, variable of interest, output file basename- PDF of ordination and top 4 ASVs, extended in the future

`bimico_summary(toy_set, "clinical_group", "ex1")`




