#' single-step preprocessing of 454 reads in BiMiCo pipeline
#'
#' Read processing function. Returns quality-filtered reads in output folder and Phyloseq object of ASVs, taxonomy and phenodata.
#' @param readfiles (Required) Path to single-end 454 demultiplexed fastq files, extension currently has to be ".fastq"
#' @param taxfile (Required) Path to taxonomy database
#' @param pheno (Required) Phenotye table
#' @param outdir (Required) Path to write quality-filtered fastq files to
#' @param trim_read_length (Optional) max. length at which to truncate reads. Default = 0 (no truncating)
#' @param mtthread (Optional) Boolean, enables multithreading (not recommended in Rstudio). Default=F
#' @keywords read processing dada2
#' @export
#' @examples
#' bimico_454()


bimico_454 <- function(readfiles, taxfile, pheno, outdir,  trim_read_length=0, mtthread=F){

prep_454(rawfqs,
         outdir,
         mtthread = mtthread,
         trim_read_length = trim_read_length)

asvtab <- asvtab_454(filtered_fqs, mtthread=mtthread)

asvtab <- asvtab[ nchar(rownames(asvtab))>=50, ]

taxtab <- asgntax(asvtab, taxfile, revcomp = T, mtthread = mtthread)

phedat <- pheno

rownames(phedat) <- gsub(".fastq", "", rownames(phedat))
colnames(asvtab) <- gsub("_filt", "", colnames(asvtab))
asvtab <- asvtab[ , order(colnames(asvtab))]
phedat <- phedat[ order(rownames(phedat)), ]

# Store data in Phyloseq object
phs <- create_phylo(asvtab, taxtab, phedat)

return(phs)

}
