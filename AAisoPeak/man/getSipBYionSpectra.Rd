% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.AAspectra.R
\name{getSipBYionSpectra}
\alias{getSipBYionSpectra}
\title{Get AAspectra object  of B and Y ions from AA sequence
with labeled SIP abundance}
\usage{
getSipBYionSpectra(AAstr, Atom = "C13", Prob = 0.01, charge = 2)
}
\arguments{
\item{AAstr}{amino acide string}

\item{Atom}{a CharacterVector C13 or N15}

\item{Prob}{C13 or N15's abundance}
}
\value{
AAspectra object
}
\description{
Get AAspectra object  of B and Y ions from AA sequence
with labeled SIP abundance
}
\examples{
getSipBYionSpectra("KHRIPCDRK", "C13", 0.05, 2)
}
