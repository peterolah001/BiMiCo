% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bimico_sumary.R
\name{bimico_summary}
\alias{bimico_summary}
\title{bimico_summary}
\usage{
bimico_summary(phylo, contrast, outfile)
}
\arguments{
\item{phylo}{(Required) Phyloseq object generated by bimico_454/bimico_illum}

\item{contrast}{(Required) column name of phenodata with the variable of interest}

\item{outfile}{(Required) output filename}
}
\description{
Generate a basic summary report based on the desired variable column of phenodata (case/control, sample source, sequencing batch...)
}
\examples{
bimico_summary()
}
\keyword{16S}
\keyword{analysis}
\keyword{data}
