#' Create ASV table from filtered PAIRED-END Illumina reads
#'
#' Given a set of PE Illumina sequencing reads pre-filtered by the 'prep_illSE' function, internally generates read error models and outputs chimera-filtered ASV table. ASV table contains ASVs in rows and samples in columns.
#' @param fwd_reads (Required) Path to PE Illumina quality-filtered fastq files
#' @param rev_reads (Required) Path to PE Illumina quality-filtered fastq files
#' @param filt_out (Required) Path to filtered reads from previous step (TrimandFilter)
#' @param rawext (Required) Raw extension, as in previous steps
#' @param mtthread (Optional) Boolean, enables multithreading (not recommended in Rstudio) Default=F
#' @param mergepairs (Optional) Boolean, to merge read pairs (recommended) or only concatenate (if there's no overlap between the majority of reads, it is recommended to inspect trimming). Default=F
#' @keywords read processing dada2
#' @export
#' @examples
#' asvtab_illPE()


asvtab_illPE <- function(fwd_reads, rev_reads, filt_out, rawext, 
                         mtthread=F,
                         mergepairs=F){
  
  
  # # store filtered fastq file names
  # ffqs <- sort(
  #   list.files(
  #     readfiles,
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
  
  
  # learn error rates: visualization of base calling error model generated for
  # sequence classification step
  
  err_F <- dada2::learnErrors(fwd_filt_fqs,
                                multithread=mtthread, nreads=100000, MAX_CONSIST = 8)
  
  err_R <- dada2::learnErrors(rev_filt_fqs,
                              multithread=mtthread, nreads=100000, MAX_CONSIST = 8)
  
  #plotErrors(err_fqs, nominalQ=TRUE)
  
  derps_F <- dada2::derepFastq(fwd_filt_fqs)
  
  derps_R <- dada2::derepFastq(rev_filt_fqs)
  
  # Apply DADA inference
  
  dada_F <- dada2::dada(derps_F,
                          err=err_F,
                          multithread=mtthread,
                          selfConsist=FALSE)
  
  dada_R <- dada2::dada(derps_R,
                        err=err_R,
                        multithread=mtthread,
                        selfConsist=FALSE)
  

  # Merge read pairs

  mergers <- dada2::mergePairs(dada_F, derps_F, 
                               dada_R, derps_R,
                               verbose=TRUE,
                               justConcatenate = mergepairs)

  
  # Create ASV table
  
  asvs <- dada2::makeSequenceTable(mergers)
  
  # remove chimeric sequences
  
  asvs.nochim <- dada2::removeBimeraDenovo(asvs,
                                           method="consensus",
                                           multithread=mtthread,
                                           verbose=TRUE)
  
  
  return(t(asvs.nochim))
  
  
}
