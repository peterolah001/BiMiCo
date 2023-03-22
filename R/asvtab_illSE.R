#' Create ASV table from filtered SINGLE-END Illumina reads
#'
#' Given a set of SE Illumina sequencing reads pre-filtered by the 'prep_illSE' function, internally generates read error models and outputs chimera-filtered ASV table. ASV table contains ASVs in rows and samples in columns.
#' @param readfiles (Required) Path to SE Illumina quality-filtered fastq DIRECTORY (e.g. filt_fqs_dir in tutorial)
#' @param n (Optional) The maximum number of reads to dereplicate at any one time. Controls peak memory requirement. Default=1e6
#' @param mtthread (Optional) Boolean, enables multithreading (not recommended in Rstudio) Default=F
#' @keywords read processing dada2
#' @export
#' @examples
#' asvtab_illSE()


asvtab_illSE <- function(readfiles, n, mtthread){
  
  
  # store filtered fastq file names
  ffqs <- sort(
    list.files(
      readfiles,
      pattern = "_filt.fastq.gz",
      full.names = TRUE
    )
  )
  
  # extract sample names
  samps <- sapply(
    strsplit(
      basename(ffqs),
      "_filt.fastq.gz"
    ),
    `[`,
    1
  )
  
  
  # learn error rates: visualization of base calling error model generated for
  # sequence classification step
  
  err_fqs <- dada2::learnErrors(ffqs,
                                multithread=mtthread, nreads=n, MAX_CONSIST = 8)
  
  #plotErrors(err_fqs, nominalQ=TRUE)
  
  # Merge identical fastq sequences
  
  derp_fqs <- dada2::derepFastq(ffqs,
                                verbose=TRUE)
  
  names(derp_fqs) <- samps
  
  # Apply DADA inference
  
  dada_fqs <- dada2::dada(derp_fqs,
                          err=err_fqs,
                          multithread=mtthread,
                          selfConsist=FALSE)
  
  # Create ASV table
  
  asvs <- dada2::makeSequenceTable(dada_fqs)
  
  # remove chimeric sequences
  
  asvs.nochim <- dada2::removeBimeraDenovo(asvs,
                                           method="consensus",
                                           multithread=FALSE,
                                           verbose=TRUE)
  
  
  return(t(asvs.nochim))
  
  
}
