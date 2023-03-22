#' standard preprocessor for PAIRED-END Illumina reads in BiMiCo pipeline
#'
#' Convenience wrapper for DADA2 read processing functions. Given a directory of raw Illumina fastq files, performs platform-specific trimming and QC, writes reads to output directory.
#'
#' @param fwd_reads (Required) Path to READ1 paired-end demultiplexed fastq files
#' @param rev_reads (Required) Path to READ1 paired-end demultiplexed fastq files
#' @param rawext (Required) READ1 fastq file precise extension, UNIQUE PART WITH READ NUMBER as string, e.g. "_1.fastq"
#' @param filt_out (Required) String specifying path to write quality-filtered fastq files to (relative or absolute path).
#' @param mtthread (Optional) Boolean, enables multithreading (not recommended in Rstudio). Default=F
#' @param trim_read_length_FWD (Optional) max. length at which to truncate reads. Default = 0 (no truncating)
#' @param trim_read_length_REV (Optional) max. length at which to truncate reads. Default = 0 (no truncating)
#' @param trim5end (Optional) trim n reads from the start of reads (unidentified adapters, barcodes). Default=0 (no trimming)
#' @keywords read processing dada2
#' @export
#' @examples
#' prep_illPE()

prep_illPE <- function(fwd_reads, rev_reads, 
                       rawext,
                       filt_out,
                       mtthread=F, 
                       trim_read_length_FWD=0,
                       trim_read_length_REV=0,
                       trim5end=0){
  
  # # store fastq file names
  # fqsF <- sort(
  #   list.files(
  #     fwd_reads,
  #     full.names = TRUE
  #   )
  # )
  # 
  # fqsR <- sort(
  #   list.files(
  #     rev_reads,
  #     full.names = TRUE
  #   )
  # )
  
  # extract sample names
  samps <- sapply(
    strsplit(
      basename(fwd_reads),
      rawext),
    `[`,
    1
  )
  
  # create dir for filtered, trimmed fastq files
  fwd_filt_fqs <- file.path(
    filt_out, "DADA2_filt_FWD",
    paste0(samps, "_df_R1.fastq.gz")
  )
  
  rev_filt_fqs <- file.path(
    filt_out, "DADA2_filt_REV",
    paste0(samps, "_df_R2.fastq.gz")
  )
  
  # default filter parameters to work with DADA2 (Illumina)
  fout <- dada2::filterAndTrim(
    fwd=fwd_reads,
    rev=rev_reads,
    filt=fwd_filt_fqs,
    filt.rev=rev_filt_fqs,
    maxN=0,
    maxEE=c(2,2),
    truncQ=2,
    truncLen = c(trim_read_length_FWD, trim_read_length_REV),
    trimLeft = trim5end,
    rm.phix=TRUE,
    compress=TRUE,
    multithread=mtthread
  )
  # Output filter summary
  
  closeAllConnections()
  
  return(fout)
  
}
