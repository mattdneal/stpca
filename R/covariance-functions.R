
#' Squared exponential covariance function
#' @export
#' @include util.R
#' @import Matrix
cov.SE <- function(X, X2, beta, D=NA, ...) {
  if (all(is.na(D))) {
    D = distanceMatrix(X, X2)
  }
  D@x = exp(beta[1]) * exp(-(D@x^2)/(2*exp(beta[2])^2))
  Matrix::diag(D) = Matrix::diag(D)+1e-9
  return(D)
}

#' Squared exponential covariance function derivatives wrt hyperparameters
#' @export
#' @include util.R
#' @import Matrix
cov.SE.d <- function(X, X2, beta, D=NA, ...) {
  if (all(is.na(D))) {
    D = distanceMatrix(X, X2)
  }
  dK1 = cov.SE(X, X2, beta, D=D)
  dK2 = (D^2)*exp(-2*beta[2]) * dK1
  return(list(dK1, dK2))
}

#' Rational quadratic covariance function
#' @export
#' @include util.R
#' @import Matrix
#' @import numDeriv
#' @examples
#' locs = matrix(c(0, 1.5), ncol=1)
#' beta = c(2, 0.5, 1.5)
#'
#' # Numerical gradients
#' grad.num  = numDeriv::grad(function(beta_) {
#'   cov.RQ(locs, beta_)[1,2]
#' }, beta)
#'
#' # Analytic gradients
#' grads = cov.RQ.d(locs, beta)
#' grad.real = vapply(1:3, function(i) grads[[i]][1,2], numeric(1))
#'
#' # Results are the same
#' stopifnot(all.equal(grad.num, grad.real))
cov.RQ <- function(X, beta, D=NA, ...) {
  # beta[1] log(sigma^2)
  # beta[2] log(lengthScale)
  # beta[3] log(alpha)
  stopifnot(length(beta)==3)
  if (all(is.na(D))) { D = distanceMatrix(X) }
  #D@x = exp(beta[1])*(1+(D@x^2)/(2*(exp(beta[2])^2)*exp(beta[3])))^(-exp(beta[3]))
  D@x = exp(beta[1])*(1 + (D@x*D@x)*0.5*exp(-2*beta[2]-beta[3]))^(-exp(beta[3]))
  return(D)
}

#' Rational quadratic covariance function derivatives wrt hyperparameters
#' @export
#' @include util.R
#' @import Matrix
cov.RQ.d <- function(X, beta, D=NA, ...) {
  stopifnot(length(beta)==3)
  if (all(is.na(D))) { D = distanceMatrix(X) }

  # Obtained using R 'grad' function
  .expr1  <- exp(beta[1])
  .expr3  <- D@x * D@x * 0.5
  .expr7  <- exp(-2 * beta[2] - beta[3])
  .expr8  <- .expr3 * .expr7
  .expr9  <- 1 + .expr8
  .expr10 <- exp(beta[3])
  .expr11 <- -.expr10
  .expr12 <- .expr9^.expr11
  .expr13 <- .expr1 * .expr12
  .expr15 <- .expr9^(.expr11 - 1)

  dK1 <- D
  dK2 <- D
  dK3 <- D
  dK1@x <- .expr13
  dK2@x <- -(.expr1 * (.expr15 * (.expr11 * (.expr3 * (.expr7 * 2)))))
  dK3@x <- -(.expr1 * (.expr12 * (log(.expr9) *
           .expr10) + .expr15 * (.expr11 * .expr8)))
  return(list(dK1, dK2, dK3))
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
