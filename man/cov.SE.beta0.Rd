% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covariance-functions.R
\name{cov.SE.beta0}
\alias{cov.SE.beta0}
\title{Computationally cheap estimate for beta0 for cov.SE.}
\usage{
cov.SE.beta0(X, locations, k)
}
\arguments{
\item{X}{The dataset being analysed with stpca}

\item{locations}{Matrix containing the location of each feature in rows}

\item{k}{Latent dimensionality used in stpca}
}
\description{
Computationally cheap estimate for beta0 for cov.SE.
}
\examples{
library(functional)
n = 10; k = 4; dim=c(10, 10); kern=Curry(cov.SE, beta=log(c(2, 0.4)))
synth = synthesize_data_kern(n, k, dim, kern, noisesd=0.2)
beta0 = cov.SE.beta0(synth$X, synth$grid, k)
stopifnot(length(beta0) == 2)
stopifnot(all(is.finite(beta0)))
}
