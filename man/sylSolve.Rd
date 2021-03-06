% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{sylSolve}
\alias{sylSolve}
\title{Solve the sylvester equation AW + WB = C for W.}
\usage{
sylSolve(A, B, C)
}
\arguments{
\item{A}{d x d positive definite matrix}

\item{B}{k x k positive definite matrix}

\item{C}{d x k matrix}
}
\value{
The solution W
}
\description{
Solve the sylvester equation AW + WB = C for W.
}
