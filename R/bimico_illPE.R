#' single-step preprocessing of paired-end Illumina reads in BiMiCo pipeline
#'
#' Read processing function. Returns quality-filtered reads in output folder and Phyloseq object of ASVs, taxonomy and phenodata.
#' @import phyloseq
#' @param fwd_reads (Required) Path to READ1 paired-end demultiplexed fastq files
#' @param rev_reads (Required) Path to READ1 paired-end demultiplexed fastq files
#' @param rawext (Required) fastq file precise extension as string, e.g. ".fastq"
#' @param taxfile (Required) Path to taxonomy database
#' @param pheno (Required) Phenotye table
#' @param filt_out (Required) Path to write quality-filtered fastq files to (subdirectories FWD and REV will be created)
#' @param outdir (Required) Path to write results to
#' @param mergepairs (Optional) Boolean, to merge read pairs (recommended) or only concatenate (if there's no overlap between the majority of reads, it is recommended to inspect trimming). Default=F
#' @param trim_read_length (Optional) max. length at which to truncate FWD and REV reads. Default = 0 (no truncating)
#' @param mtthread (Optional) Boolean, enables multithreading (not recommended in Rstudio). Default=F
#' @keywords read processing dada2
#' @export
#' @examples
#' bimico_illPE()


bimico_illPE <- function(fwd_reads, rev_reads, rawext, 
                         taxfile, pheno, filt_out, outdir, mergepairs=F,
                         trim_read_length=0, mtthread=F){

  dir.create(outdir)
  
  prep_illPE(fwd_reads=fwd_reads, 
                         rev_reads=rev_reads, 
                         rawext=rawext,
                         filt_out=filt_out,
                         mtthread=F, 
                         trim_read_length_FWD=trim_read_length,
                         trim_read_length_REV=trim_read_length,
                         trim5end=0)

  
asvtab <- asvtab_illPE(fwd_reads = fwd_reads,
                       rev_reads = rev_reads,
                       filt_out = filt_out,
                       rawext = rawext,
                       mtthread=mtthread,
                       mergepairs = mergepairs)

# asvtab <- asvtab[ nchar(rownames(asvtab))>=50, ]

taxtab <- asgntax(asvtab, taxfile, revcomp = T, mtthread = mtthread)

phedat <- pheno

 colnames(asvtab) <- gsub("_df_R1.fastq.gz", "", colnames(asvtab))
asvtab <- asvtab[ , order(colnames(asvtab))]
phedat <- phedat[ order(rownames(phedat)), ]

print("samples matching between phenotable and ASV table:    ") 
print(summary(rownames(phedat) %in% colnames(asvtab)))

saveRDS(asvtab, file = file.path(outdir,"asv_table.RDS"))
saveRDS(taxtab, file = file.path(outdir,"tax_table.RDS"))

# Store data in Phyloseq object
phs <- create_phylo(asvtab, taxtab, phedat)

return(phs)

}
