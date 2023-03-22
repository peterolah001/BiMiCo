## TODO: define optional args
#       reads file name not strictly "fastq"!
#       make chimera removal optional

#' Create ASV table from filtered 454 reads
#'
#' Given a set of 454 sequencing reads pre-filtered by the 'prep_454' function, internally generates read error models and outputs chimera-filtered ASV table. ASV table contains ASVs in rows and samples in columns.
#' @param readfiles (Required) Path to 454 quality-filtered fastq files
#' @param n (Optional) The maximum number of reads to dereplicate at any one time. Controls peak memory requirement. Default=1e6
#' @param mtthread (Optional) Boolean, enables multithreading (not recommended in Rstudio) (default=F)
#' @param pattern (Optional) Pattern (e.g., extension) in input fastq filename, only files with this pattern will be read (default="_filt.fastq.gz")
#' @param extname (Optional) Suffix of input files. (default="_filt.fastq.gz")
#' @keywords read processing dada2
#' @export
#' @examples
#' asvtab_454()


asvtab_454 <- function(readfiles, n, mtthread = F, pattern = "_filt.fastq.gz", extname="_filt.fastq.gz"){

  # store filtered fastq file names
  ffqs <- sort(
    list.files(
      readfiles,
      pattern = pattern,
      full.names = TRUE
    )
  )
  
  # extract sample names
  samps <- sapply(
    strsplit(
      basename(ffqs),
      extname
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
                        HOMOPOLYMER_GAP_PENALTY=-1,
                        BAND_SIZE=32,
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
