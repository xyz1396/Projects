% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.AAspectra.R
\name{getPrecursorSpectra}
\alias{getPrecursorSpectra}
\title{Get AAspectra object of precursor from AA sequence with natural SIP abundance}
\usage{
getPrecursorSpectra(AAstr, charge = 1)
}
\arguments{
\item{AAstr}{amino acide string}

\item{charge}{charge of the ion}
}
\value{
AAspectra object
}
\description{
Get AAspectra object of precursor from AA sequence with natural SIP abundance
}
\examples{
getPrecursorSpectra("KHRIP",2)
}
