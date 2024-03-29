% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep_illPE.R
\name{prep_illPE}
\alias{prep_illPE}
\title{standard preprocessor for PAIRED-END Illumina reads in BiMiCo pipeline}
\usage{
prep_illPE(fwd_reads, rev_reads, rawext, filt_out, mtthread = F,
  trim_read_length_FWD = 0, trim_read_length_REV = 0, trim5end = 0)
}
\arguments{
\item{fwd_reads}{(Required) Path to READ1 paired-end demultiplexed fastq files}

\item{rev_reads}{(Required) Path to READ1 paired-end demultiplexed fastq files}

\item{rawext}{(Required) READ1 fastq file precise extension, UNIQUE PART WITH READ NUMBER as string, e.g. "_1.fastq"}

\item{filt_out}{(Required) String specifying path to write quality-filtered fastq files to (relative or absolute path).}

\item{mtthread}{(Optional) Boolean, enables multithreading (not recommended in Rstudio). Default=F}

\item{trim_read_length_FWD}{(Optional) max. length at which to truncate reads. Default = 0 (no truncating)}

\item{trim_read_length_REV}{(Optional) max. length at which to truncate reads. Default = 0 (no truncating)}

\item{trim5end}{(Optional) trim n reads from the start of reads (unidentified adapters, barcodes). Default=0 (no trimming)}
}
\description{
Convenience wrapper for DADA2 read processing functions. Given a directory of raw Illumina fastq files, performs platform-specific trimming and QC, writes reads to output directory.
}
\examples{
prep_illPE()
}
\keyword{dada2}
\keyword{processing}
\keyword{read}
