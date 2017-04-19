#' @include util.R

#' Squared exponential covariance function
#' @export
cov.SE <- function(X, X2, beta, D=NA, ...) {
  if (all(is.na(D))) {
    D = distanceMatrix(X, X2)
  }
  D@x = exp(beta[1]) * exp(-(D@x^2)/(2*exp(beta[2])^2))
  Matrix::diag(D) = Matrix::diag(D)+1e-9
  return(D)
}

#' Squared exponential covariance function partial derivatives wrt hyperparameters
#' @export
cov.SE.d <- function(X, X2, beta, D=NA, ...) {
  if (all(is.na(D))) {
    D = distanceMatrix(X, X2)
  }
  dK1 = cov.SE(X, X2, beta)
  dK2 = (D^2)*exp(-2*beta[2]) * dK1
  return(list(dK1, dK2))
}

#' Rational quadratic covariance function
#' @export
cov.RQ <- function(X, beta, ...) {
  # beta[1] log(sigma^2)
  # beta[2] log(lengthScale)
  # beta[3] log(alpha)
  D   = distanceMatrix(X)
  D@x = exp(beta[1])*(1+(D@x^2)/(2*(exp(beta[2])^2)*exp(beta[3])))^(-exp(beta[3]))
  return(D)
}

cov.independent <- function(X, X2=NA, beta=c(), D=NA, max.dist=NA) {
  stopifnot(nrow(X) == length(beta))
  return(sparseMatrix(i=1:nrow(X), j=1:nrow(X), x=exp(beta), dims=c(nrow(X), nrow(X))))
}

cov.independent.d <- function(X, X2=NA, beta=c(), D=NA, max.dist=NA) {
  return(lapply(1:nrow(X), function(i) {
    sparseMatrix(i=i, j=i, x=1, dims=c(nrow(X), nrow(X)))
  }))
}

cov.triangular <- function(X, beta, ...) {
  D   = distanceMatrix(X, max.dist=exp(beta[2]))
  D@x = exp(beta[1])*(1 - D@x/exp(beta[2]))
  return(D)
}

cov.triangular.d <- function(X, beta, ...) {
  D   = distanceMatrix(X, max.dist=exp(beta[2]))

  dKdBeta1 = D
  dKdBeta1@x = 1 - D@x/exp(beta[2])

  dKdBeta2 = D
  dKdBeta2@x = exp(beta[1])*D@x/(exp(beta[2])^2)
  return(list(dKdBeta1, dKdBeta2))
}
