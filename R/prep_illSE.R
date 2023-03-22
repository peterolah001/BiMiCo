#' standard preprocessor for SINGLE-END Illumina reads in BiMiCo pipeline
#'
#' Convenience wrapper for DADA2 read processing functions. Given a directory of raw Illumina single-end fastq files, performs platform-specific trimming and QC, writes reads to output directory.
#'
#' @param readfiles (Required) Path to single-end Illumina demultiplexed fastq files
#' @param outdir (Required) String specifying output directory for processed files; will be created as subdirectory of input dir.
#' @param file_suffix (Optional) Please specify if input files have special suffixes other than ".fastq", e.g. "_filt.hf.fq". Default=".fastq"
#' @param mtthread (Optional) Boolean, enables multithreading (not recommended in Rstudio). Default=F
#' @param trim_read_length (Optional) max. length at which to truncate reads. Default = 0 (no truncating)
#' @param trim5end (Optional) trim n reads from the start of reads (unidentified adapters, barcodes). Default=0 (no trimming)
#' @keywords read processing dada2
#' @export
#' @examples
#' prep_illSE()

prep_illSE <- function(readfiles, outdir, file_suffix=".fastq", mtthread=F, trim_read_length=0, trim5end=0){

  # # store fastq file names
  # fqs <- sort(
  #   list.files(
  #     readfiles,
  #     full.names = TRUE
  #     )
  #   )

  # extract sample names
  samps <- sapply(
    strsplit(
    basename(readfiles),
    file_suffix),
    `[`,
    1
    )

  # outdir
  dir.create(outdir, recursive = T)
  
  # create paths for filtered, trimmed fastq files
  filt_fqs <- file.path(
    outdir,
    paste0(samps, "_filt.fastq.gz")
  )


  # default filter parameters to work with DADA2 (Illumina)
  fout <- dada2::filterAndTrim(
    readfiles,
    filt_fqs,
                        maxN=0,
                        maxEE=2,
                        truncQ=2,
                        truncLen = trim_read_length,
                        trimLeft = trim5end,
                        rm.phix=TRUE,
                        compress=TRUE,
                        multithread=mtthread
    )
# Output filter summary

  closeAllConnections()

  return(fout)

  }
