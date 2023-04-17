#' Create ASV table from filtered PAIRED-END Illumina reads
#'
#' Given a set of PE Illumina sequencing reads pre-filtered by the 'prep_illSE' function, internally generates read error models and outputs chimera-filtered ASV table. ASV table contains ASVs in rows and samples in columns.
#' @param filt_out (Required) Path to filtered reads from previous step (dir containing "DADA2_filt_FWD" and DADA2_filt_REV" dirs)
#' @param batch (Required) List indicating the sequencing batch (e.g., a column of phenodata with batch information)
#' @param mtthread (Optional) Boolean, enables multithreading (not recommended in Rstudio) Default=F
#' @param mergepairs (Optional) Boolean, to merge read pairs (recommended) or only concatenate (if there's no overlap between the majority of reads, it is recommended to inspect trimming). Default=F
#' @keywords read processing dada2
#' @export
#' @examples
#' asvtab_illPE()

# TODO: make batch and output destinations optional

asvtab_illPE <- function(filt_out, batch,
                         mtthread=F,
                         mergepairs=F){
  
  asvnochims <- list()
  
  for (i in levels(as.factor(batch))) {
  
  # Get fwd and rev filtered fqs
    fwd_filt_fqs <- sort(list.files(filt_out,pattern="_df_R1.fastq.gz", recursive = T, full.names = T))
    rev_filt_fqs <- sort(list.files(filt_out,pattern="_df_R2.fastq.gz", recursive = T, full.names = T))
  
  # get files of current batch
    fwd_filt_fqs <- cbind(as.data.frame(fwd_filt_fqs),batch)
    fwd_filt_fqs <- as.character(fwd_filt_fqs[ fwd_filt_fqs$batch==i, 1])
  
    rev_filt_fqs <- cbind(as.data.frame(rev_filt_fqs),batch)
    rev_filt_fqs <- as.character(rev_filt_fqs[ rev_filt_fqs$batch==i, 1])
    
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
  
  asvnochims[[i]] <-  as.data.frame(asvs.nochim) 
  asvcols <- unlist(lapply(asvnochims, rownames))
  
  }
  
  asv_table <-  Reduce(function (...) { merge(..., all = TRUE) },   # Full join
                       asvnochims) 
   rownames(asv_table) <- asvcols
  asv_table[ is.na(asv_table) ] <- 0
  
  asv_table <- t(asv_table)

return(asv_table)
#return(t(asvs.nochim))
  
}
