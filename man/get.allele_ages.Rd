% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/manipulate_sele.R
\name{get.allele_ages}
\alias{get.allele_ages}
\title{Compile allele ages file}
\usage{
get.allele_ages(mut, freq, sele)
}
\arguments{
\item{mut}{data.table. mut file}

\item{freq}{data.frame. freq file}

\item{sele}{data.frame. sele file}
}
\value{
Returns a data table.
}
\description{
This function aggregates mut, freq, and sele files inferred by Relate into an allele_ages file.
}
\examples{
mut         <- read.mut(system.file("extdata/example.mut.gz", package = "relater"), CHR = 1)
sele        <- read.sele(system.file("extdata/example.sele.gz", package = "relater"), CHR = 1)
freq        <- read.freq(system.file("extdata/example.freq.gz", package = "relater"), CHR = 1)
get.allele_ages(mut, freq, sele)
}
