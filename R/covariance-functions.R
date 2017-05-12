#' Squared exponential covariance function
#' @export
#' @include util.R
#' @import Matrix
#' @examples
#' grid = matrix(1:10, ncol=1)
#' beta = rnorm(2)
#' K    = cov.SE(grid, beta=beta)
#' stopifnot(all(Matrix::diag(K)==exp(beta[1])))
cov.SE <- function(X, X2, beta, D=NA, ...) {
  if (all(is.na(D))) {
    D = distanceMatrix(X, X2)
  }
  #D@x = exp(beta[1]) * exp(-0.5*(D@x^2)*exp(-2*beta[2]))
  D = exp(beta[1]) * exp(-0.5*(D^2)*exp(-2*beta[2]))
  #Matrix::diag(D) = Matrix::diag(D)+1e-12
  return(D)
}

#' Squared exponential covariance function derivatives wrt hyperparameters
#' @export
#' @include util.R
#' @import Matrix
#' @import numDeriv
#' @examples
#' point1 = matrix(rnorm(1), ncol=1)
#' point2 = matrix(rnorm(1), ncol=1)
#' beta   = rnorm(2) # Logarithms of variance and length scale
#'
#' Ks.1point = cov.SE.d(point1, beta=beta)
#'
#' # Derivative wrt variance at zero distance should always be exp(beta[1])
#' stopifnot(all.equal(Ks.1point[[1]][1,1], exp(beta[1])))
#'
#' # Derivative wrt lengthscale at zero distance should always be 0
#' stopifnot(all.equal(Ks.1point[[2]][1,1], 0))
#'
#' # Identical tests with numerical gradient
#' library(numDeriv)
#' Ks.1point.num = grad(function(beta_) {
#'   cov.SE(point1, beta=beta_)[1,1]
#' }, x=beta)
#' stopifnot(all.equal(Ks.1point.num[1], exp(beta[1])))
#' stopifnot(all.equal(Ks.1point.num[2], 0))
#'
#' Ks.2points = cov.SE.d(point1, point2, beta=beta)
#' Ks.2points.num = grad(function(beta_) {
#'   as.numeric(cov.SE(point1, point2, beta=beta_))
#' }, x=beta)
#'
#' # Check numerical gradient equals analytic gradient
#' stopifnot(all.equal(as.numeric(Ks.2points[[1]]), Ks.2points.num[1]))
#' stopifnot(all.equal(as.numeric(Ks.2points[[2]]), Ks.2points.num[2]))
cov.SE.d <- function(X, X2, beta, D=NA, ...) {
  if (all(is.na(D))) {
    D = distanceMatrix(X, X2)
  }
  dK1 = exp(beta[1])*exp(-0.5*(D^2)*exp(-2*beta[2]))
  dK2 = exp(beta[2])*(D^2)*exp(-3*beta[2]) * dK1
  return(list(dK1, dK2))
}

#' Computationally cheap estimate for beta0 for cov.SE.
#' @export
#' @include synthesize-data.R
#' n = 10; k = 4; dim=c(10, 10); kern=Curry(cov.SE, beta=log(c(2, 0.4)))
#' synth = synthesize_data_kern(n, k, dim, kern, noisesd=0.2)
#' beta0 = cov.SE.beta0(synth$X, synth$grid)
#' stopifnot(all(is.finite(beta0)))
cov.SE.beta0 <- function(X, locations) {
  sigSqf0 = mean(apply(X, 2, var)/k)
  Rsq    = apply(locations, 1, function(loc) colSums((t(locations)-loc)^2))
  C      = cov(as.matrix(X))

  # Upper-bounded by the maximum distance, since we cant learn a much larger distance than this!
  l0     = min(mean(sqrt(0.5*Rsq/(log(k*sigSqf0) - log(mean(C)))), na.rm=TRUE),
               sqrt(max(Rsq)))
  beta0 = log(c(sigSqf0, l0))

  return(beta0)
}

#' Rational quadratic covariance function
#' @export
#' @include util.R
#' @import Matrix
#' @import numDeriv
#' @examples
#' locations = matrix(rnorm(10), ncol=2)
#' beta      = rnorm(3) # Logarithms of variance, length scale & alpha
#'
#' K = cov.RQ(locations, beta=beta)
#' stopifnot(all(Matrix::diag(K)==exp(beta[1]))) # Diagonal is exp(beta[1])
#' stopifnot(all(svd(K)$d>0)) # K is positive definite
#' stopifnot(all(K[upper.tri(K)]<K[1,1])) # Largest element is on diagonal
cov.RQ <- function(X, X2, beta, D=NA, ...) {
  stopifnot(length(beta)==3)
  if (all(is.na(D))) { D = distanceMatrix(X, X2) }
  D = exp(beta[1])*(1 + (D*D)*0.5*exp(-2*beta[2]-beta[3]))^(-exp(beta[3]))
  return(D)
}

#' Rational quadratic covariance function derivatives wrt hyperparameters
#' @export
#' @include util.R
#' @import Matrix
#' @import numDeriv
#' @examples
#' point1 = matrix(rnorm(1), ncol=1)
#' point2 = matrix(rnorm(1), ncol=1)
#' beta   = rnorm(3) # Logarithms of variance, length scale & alpha
#'
#' library(numDeriv)
#' Ks.2points = vapply(cov.RQ.d(point1, point2, beta=beta), function(K) {
#'   K[1,1]
#' }, numeric(1))
#' Ks.2points.num = grad(function(beta_) {
#'   as.numeric(cov.RQ(point1, point2, beta=beta_))
#' }, x=beta)
#' stopifnot(all.equal(Ks.2points, Ks.2points.num))
cov.RQ.d <- function(X, X2, beta, D=NA, ...) {
  stopifnot(length(beta)==3)
  if (all(is.na(D))) { D = distanceMatrix(X, X2) }

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

#' Computationally cheap estimate for beta0 for cov.RQ.
#' @export
#' @include synthesize-data.R
#' n = 10; k = 4; dim=c(10, 10); kern=Curry(cov.SE, beta=log(c(2, 0.4, 0.3)))
#' synth = synthesize_data_kern(n, k, dim, kern, noisesd=0.2)
#' beta0 = cov.RQ.beta0(synth$X, synth$grid)
#' stopifnot(all(is.finite(beta0)))
cov.RQ.beta0 <- function(X, locations) {
  sigSqf0 = mean(apply(X, 2, var)/k)
  Rsq    = apply(locations, 1, function(loc) colSums((t(locations)-loc)^2))
  C      = cov(as.matrix(X))

  # Upper-bounded by the maximum distance, since we cant learn a much larger distance than this!
  l0     = min(mean(sqrt(0.5*Rsq/(log(k*sigSqf0) - log(mean(C)))), na.rm=TRUE),
               sqrt(max(Rsq)))

  #browser()
  #f = function(alpha_) {
  #  (1 + Rsq/(2*alpha_*l0*l0))^(-alpha_) - C
  #}
  #alpha0 = uniroot(f, interval=c(0, 10))
  #alpha0 = uniroot(f, interval=c(0.01, 2))$root
  alpha0 = 2

  beta0 = log(c(sigSqf0, l0, alpha0))

  return(beta0)
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
