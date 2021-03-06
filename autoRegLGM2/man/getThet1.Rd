% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getThet1.R
\name{getThet1}
\alias{getThet1}
\title{Obtaining Initial trait estimates}
\usage{
getThet1(zMat, alpha, ak, bk, varL1)
}
\arguments{
\item{zMat}{Latent Normally distributed outcome variables}

\item{alpha}{initial ability estimates from the HLM}

\item{ak}{A numeric vector of item-level discrimination values}

\item{bk}{A numeric vector of item-level difficulties}

\item{nItems}{Number of items}
}
\description{
Before an FFBS is run, we need an initial estimate of the trait.
This function uses values from getZs(), the HLM, and item characteristics to get this initial estimate.
It returns a list of precisions and means for each person.
}
\examples{
getThet1()
}
\keyword{Initial}
\keyword{initializing}

