% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evidence-approximation.R
\name{stpca.log_prior_d}
\alias{stpca.log_prior_d}
\title{Compute the partial derivatives of the log prior with respect to the
hyperparameters with the value of beta provided.}
\usage{
stpca.log_prior_d(W, beta, K, dK)
}
\arguments{
\item{W}{Loadings matrix}

\item{beta}{hyperparameter values}

\item{K}{Prior covariance matrix}

\item{dK}{Prior covariance matrix derivatives}
}
\value{
Partial derivatives of log prior
}
\description{
Compute the partial derivatives of the log prior with respect to the
hyperparameters with the value of beta provided.
}
