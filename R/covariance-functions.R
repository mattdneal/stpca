#' Noisy RQ covariance function
#'
#' The sum of the Rational Quadratic and the independent covariance functions.
#' Interpreted as the smooth RQ covariance function plus iid noise. This is
#' implemented as the sum of two covariance functions, with the first three
#' elements of beta being sent to the RQ covariance function, and the last
#' element being sent to the independent covariance function.
#'
#' @param X Matrix of data
#' @param X2 (optional) second matrix of data; if omitted, X is used.
#' @param beta Hyperparameters; beta[1] is the log signal variance, beta[2] is
#'          the log length scale, beta[3] is log alpha, and beta[4] is the
#'          variance of the noise.
#' @export
#' @include util.R
#' @import Matrix
cov.noisy.RQ <- function(X, X2, beta, D=NA, ...) {
  stopifnot(is.matrix(X))
  stopifnot(length(beta)==4)
  if (all(is.na(D))) { D = distanceMatrix(X, X2) }
  K = cov.RQ(X, X2, beta[1:3], D) + cov.independent(X, X2, beta[4], D)
  return(K)
}

#' Partial derivatives of the noisy RQ covariance function
#'
#' @param X Matrix of data
#' @param X2 (optional) second matrix of data; if omitted, X is used.
#' @param beta Hyperparameters; beta[1] is the log signal variance, beta[2] is
#'          the log length scale, beta[3] is log alpha, and beta[4] is the
#'          variance of the noise.
#' @export
#' @include util.R
#' @import Matrix
cov.noisy.RQ.d <- function(X, X2, beta, D=NA, ...) {
  stopifnot(is.matrix(X))
  stopifnot(length(beta)==4)
  if (all(is.na(D))) { D = distanceMatrix(X, X2) }
  return(append(cov.RQ.d(X, X2, beta[1:3], D),
                cov.independent.d(X, X2, beta[4], D)))
}

#' Noisy SE covariance function
#'
#' The sum of the Squared Exponential and the independent covariance functions.
#' Interpreted as the smooth SE covariance function plus iid noise. This is
#' implemented as the sum of two covariance functions, with the first two
#' elements of beta being sent to the SE covariance function, and the last
#' element being sent to the independent covariance function.
#'
#' @param X Matrix of data
#' @param X2 (optional) second matrix of data; if omitted, X is used.
#' @param beta Hyperparameters; beta[1] is the log signal variance, beta[2] is the log length scale, beta[3] is the variance of the noise.
#' @export
#' @include util.R
#' @import Matrix
cov.noisy.SE <- function(X, X2, beta, D=NA, ...) {
  stopifnot(is.matrix(X))
  stopifnot(length(beta)==3)
  if (all(is.na(D))) {
    D = distanceMatrix(X, X2)
  }
  K = (cov.SE(X, X2, beta=beta[1:2], D=D) +
       cov.independent(X, X2, beta=beta[3], D=D))
  return(K)
}

#' Partial derivatives of the noisy SE covariance function
#'
#' @param X Matrix of data
#' @param X2 (optional) second matrix of data; if omitted, X is used.
#' @param beta Hyperparameters; beta[1] is the log signal variance, beta[2] is the log length scale, beta[3] is the variance of the noise.
#' @export
#' @include util.R
#' @import Matrix
cov.noisy.SE.d <- function(X, X2, beta, D=NA, ...) {
  stopifnot(is.matrix(X))
  stopifnot(length(beta)==3)
  if (all(is.na(D))) { D = distanceMatrix(X, X2) }
  return(append(cov.SE.d(X, X2, beta[1:2], D),
                cov.independent.d(X, X2, beta[3], D)))
}

#' Computationally cheap estimate for beta0 for cov.noisy.SE.
#' @param X The dataset being analysed with stpca
#' @param locations Matrix containing the location of each feature in rows
#' @param k Latent dimensionality used in stpca
#' @export
#' @include synthesize-data.R
cov.noisy.SE.beta0 <- function(X, locations, k) {
  stopifnot(is.matrix(locations))
  n = nrow(X); d = ncol(X)

  covar.svd = svd(scale(X, scale=FALSE)/sqrt(n), nu=0, nv=0)
  covar.eigval = covar.svd$d^2
  sigSq = sum(covar.eigval[-(1:k)])/(d-k)

  # \sigma^2_f <- 1/k mean( diag(cov(X) - \sigma^2\mathit{I}) )
  sigSqk0 = mean((apply(X, 2, var) - sigSq)/k)
  Rsq    = apply(locations, 1, function(loc) colSums((t(locations)-loc)^2))
  C      = cov(as.matrix(X))

  # Upper-bounded by the maximum distance, since we cant learn a much larger distance than this!
  l0    = mean(sqrt(0.5*Rsq/(log(k*sigSqk0) - log(mean(C)))), na.rm=TRUE)

  # Can't learn lengthscales much longer than longest or shorter than shortest
  # observed distance.
  l0    = min(l0, 0.9*sqrt(max(Rsq)))
  l0    = max(l0, 1.1*sqrt(min(Rsq)))
  beta0 = log(c("logsigSqk"=sigSqk0, "logl"=l0, "logsigSqn"=sigSq))

  return(beta0)
}

#' Squared exponential covariance function.
#'
#' The squared exponential covariance function. This produces a semidefinite
#' covariance matrix, and should only be used when constructing new covariance
#' functions. E.g., the squared exponential plus independant.
#'
#' @param X Matrix of data
#' @param X2 (optional) second matrix of data; if omitted, X is used.
#' @param beta Hyperparameters; beta[1] is the log signal variance, beta[2] is the log length scale.
#' @export
#' @include util.R
#' @import Matrix
#' @examples
#' # Confirm that diagonal of covariance matrix is the specified variance, and
#' # the off-diagonals are less than this.
#' grid = matrix(1:10, ncol=1)
#' beta = rnorm(2)
#' K    = cov.SE(grid, beta=beta)
#' stopifnot(all(Matrix::diag(K)==exp(beta[1])))
#' stopifnot(all(K[upper.tri(K)]<exp(beta[1])))
cov.SE <- function(X, X2, beta, D=NA, ...) {
  stopifnot(is.matrix(X))
  stopifnot(length(beta)==2)
  if (all(is.na(D))) {
    D = distanceMatrix(X, X2)
  }
  D@x = exp(beta[1]) * exp(-0.5*(D@x*D@x)*exp(-2*beta[2]))
  #D = exp(beta[1]) * exp(-0.5*(D^2)*exp(-2*beta[2]))
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
  stopifnot(is.matrix(X))
  if (all(is.na(D))) {
    D = distanceMatrix(X, X2)
  }
  dK1 = D
  dK1@x = exp(beta[1])*exp(-0.5*(D@x*D@x)*exp(-2*beta[2]))
  dK2 = D
  dK2@x = exp(beta[2])*(D@x*D@x)*exp(-3*beta[2]) * dK1@x
  return(list(dK1, dK2))
}

#' Computationally cheap estimate for beta0 for cov.SE.
#' @param X The dataset being analysed with stpca
#' @param locations Matrix containing the location of each feature in rows
#' @param k Latent dimensionality used in stpca
#' @export
#' @include synthesize-data.R
#' @examples
#' library(functional)
#' n = 10; k = 4; dim=c(10, 10); kern=Curry(cov.SE, beta=log(c(2, 0.4)))
#' synth = synthesize_data_kern(n, k, dim, kern, noisesd=0.2)
#' beta0 = cov.SE.beta0(synth$X, synth$grid, k)
#' stopifnot(length(beta0) == 2)
#' stopifnot(all(is.finite(beta0)))
cov.SE.beta0 <- function(X, locations, k) {
  return(cov.noisy.SE.beta0(X, locations, k)[c("logsigSqk", "logl")])
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
  stopifnot(is.matrix(X))
  stopifnot(length(beta)==3)
  if (all(is.na(D))) { D = distanceMatrix(X, X2) }
  D@x = exp(beta[1])*(1 + (D@x*D@x)*0.5*exp(-2*beta[2]-beta[3]))^(-exp(beta[3]))
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
  stopifnot(is.matrix(X))
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
#' @param X The dataset being analysed with stpca
#' @param locations Matrix containing the location of each feature in rows
#' @param k Latent dimensionality used in stpca
#' @export
#' @include synthesize-data.R
#' @examples
#' # Construct some synthetic data, initialise beta for the RQ kernel from this dataset.
#' library(functional)
#' n = 10; k = 4; dim=c(10, 10); kern=Curry(cov.SE, beta=log(c(2, 0.4)))
#' synth = synthesize_data_kern(n, k, dim, kern, noisesd=0.2)
#' beta0 = cov.RQ.beta0(synth$X, synth$grid, k)
#' stopifnot(length(beta0) == 3)
#' stopifnot(all(is.finite(beta0)))
cov.RQ.beta0 <- function(X, locations, k) {
  stopifnot(is.matrix(locations))
  beta0SE = cov.SE.beta0(X, locations, k)

  #browser()
  #f = function(alpha_) {
  #  (1 + Rsq/(2*alpha_*l0*l0))^(-alpha_) - C
  #}
  #alpha0 = uniroot(f, interval=c(0, 10))
  #alpha0 = uniroot(f, interval=c(0.01, 2))$root
  alpha0 = 2

  beta0 = c(beta0SE, "logalpha0"=log(alpha0))

  return(beta0)
}

#' Computationally cheap estimate for beta0 for cov.noisy.RQ.
#' @param X The dataset being analysed with stpca
#' @param locations Matrix containing the location of each feature in rows
#' @param k Latent dimensionality used in stpca
#' @export
#' @include synthesize-data.R
cov.noisy.RQ.beta0 <- function(X, locations, k) {
  beta0SE = cov.noisy.SE.beta0(X, locations, k)
  beta0   = c(beta0SE, "logalpha0"=log(alpha0))
  return(beta0)
}

#' Independant covariance function. Is zero everywhere except for inputs with
#' zero distance. Has single hyperparameter of log signal variance. All nonzero
#' outputs are equal to the signal variance.
#'
#' @param X Matrix of data
#' @param X2 (optional) second matrix of data; if omitted, X is used.
#' @param beta The single hyperparameter: the log of the signal variance
#' @export
#' @include util.R
#' @import Matrix
cov.independent <- function(X, X2, beta, D=NA, ...) {
  stopifnot(is.matrix(X))
  stopifnot(length(beta)==1)
  if (all(is.na(D))) {
    D = distanceMatrix(X, X2, max.dist=1e-8)
  }
  D@x=1/D@x
  D = as(D, "TsparseMatrix")
  keep = which(is.infinite(D@x))
  return(sparseMatrix(i=D@i[keep], j=D@j[keep], x=rep(exp(beta), length(keep)),
         dims=c(nrow(D), ncol(D)), index1=FALSE, symmetric=TRUE))
}

#' Derivative of the independent covariance function. Does not depend
#' on the value of \eqn{\sigma^2_k}; is always 1 everywhere the inputs have
#' zero distance, and zero everywhere else.
#'
#' @param X Matrix of data
#' @param X2 (optional) second matrix of data; if omitted, X is used.
#' @param beta The single hyperparameter: the log of the signal variance
#' @export
#' @include util.R
#' @import Matrix
cov.independent.d <- function(X, X2, beta, D=NA, ...) {
  stopifnot(is.matrix(X))
  stopifnot(length(beta)==1)
  if (all(is.na(D))) {
    D = distanceMatrix(X, X2, max.dist=1e-8)
  }
  D@x=1/D@x
  D = as(D, "TsparseMatrix")
  keep = which(is.infinite(D@x))
  return(list(as(
    sparseMatrix(i=D@i[keep], j=D@j[keep], x=rep(exp(beta), length(keep)),
                 dims=c(nrow(D), ncol(D)), index1=FALSE, symmetric=TRUE),
    "symmetricMatrix")))
}

#' Noisy MR
#'
#' @param X Matrix of data
#' @param beta Hyperparameters; beta[1] is the log signal variance,
#'   beta[-1] are the ncol(X) length scales.
#' @export
#' @include util.R
#' @import Matrix
cov.MR <- function(X, beta, ...) {
  stopifnot(is.matrix(X))
  stopifnot(length(beta)==ncol(X)+1)
  sigVar = exp(beta[1])  # Signal variance: k(0)
  supLen = exp(beta[-1]) # Support length; if d>supLen, k(d)=0

  D = distanceMatrix(X%*%diag(1/supLen), max.dist=1)

  r = D@x
  D@x = sigVar*((2+cos(2*pi*r))*(1-r)/3 + sin(2*pi*r)/(2*pi))
  return(D)
}

#' Partial derivatives of the MR
#'
#' @param X Matrix of data
#' @param beta Hyperparameters; beta[1] is the log signal variance,
#'   beta[-1] are the ncol(X) length scales.
#' @export
#' @include util.R
#' @import Matrix
cov.MR.d <- function(X, beta, ...) {
  stopifnot(is.matrix(X))
  stopifnot(length(beta)==ncol(X)+1)
  sigVar = exp(beta[1])  # Signal variance: k(0)
  supLen = exp(beta[-1]) # Support length; if d>supLen, k(d)=0

  D = distanceMatrix(X%*%diag(1/supLen), max.dist=1)

  r = D@x

  dK1x  = ((2+cos(2*pi*r))*(1-r)/3 + sin(2*pi*r)/(2*pi))*sigVar

  dK1 = D; dK1@x = dK1x
  derivs = list(logSigVar=dK1)

  ij = getij(D)
  for (dimension in seq_len(ncol(X))) {
    X2 = X[,dimension]/supLen[dimension]
    r2 = (X2[ij$i] - X2[ij$j])^2

    dK2x = ((pi*(1-r)*cos(pi*r) + sin(pi*r)) *
            sin(pi*r)*r2/(r*supLen[dimension])) *
           4*sigVar*supLen[dimension]/3
    dK2   = D
    dK2@x = dK2x
    diag(dK2) = 0
    derivs[[paste("logl", dimension, sep='')]] = dK2
  }

  return(derivs)
}

#' Noisy MR covariance function
#'
#' @param X Matrix of data
#' @param beta Hyperparameters; beta[1] is the log signal variance,
#'   beta[2:(ncol(X)+1)] are the ncol(X) length scales. beta[ncol(X)+2]
#'   is the noise variance
#' @export
#' @include util.R
#' @import Matrix
cov.noisy.MR <- function(X, beta, ...) {
  stopifnot(is.matrix(X))
  stopifnot(length(beta)==ncol(X)+2)
  K = cov.MR(X, beta=beta[1:(ncol(X)+1)]) +
      cov.independent(X, beta=beta[-(1:(ncol(X)+1))])
  return(K)
}

#' Partial derivatives of the noisy MR covariance function
#'
#' @param X Matrix of data
#' @param beta Hyperparameters; beta[1] is the log signal variance,
#'   beta[2:(ncol(X)+1)] are the ncol(X) length scales. beta[ncol(X)+2]
#'   is the noise variance
#' @export
#' @include util.R
#' @import Matrix
cov.noisy.MR.d <- function(X, beta, ...) {
  stopifnot(is.matrix(X))
  stopifnot(length(beta)==ncol(X)+2)
  return(append(cov.MR.d(X, beta=beta[1:(ncol(X)+1)]),
                cov.independent.d(X, beta=beta[-(1:(ncol(X)+1))])))
}

#' Computationally cheap estimate for beta0 for cov.MR
#' @param X The dataset being analysed with stpca
#' @param locations Matrix containing the location of each feature in rows
#' @param k Latent dimensionality used in stpca
#' @export
cov.MR.beta0 <- function(X, locations, k) {
  stopifnot(is.matrix(locations))
  beta0 = cov.SE.beta0(X, locations, k)
  beta0[2] = beta0[2] + log(2.5)
  beta0[2:(ncol(locations)+1)] = beta0[2]
  names(beta0)[2:(ncol(locations)+1)] = paste("logSupport",
    1:ncol(locations), sep='')
  return(beta0)
}

#' Taper a covariance function
#'
#' Takes a covariance function and tapers it by returning a new covariance
#' function which is the original multiplied by the MR covariance function. The
#' support length of the new covariance function will thus be the support
#' length of the MR (assuming that the original is not compactly supported).
#'
#' @param cov The covariance function to be tapered
#' @param supLen The support length of the MR
#' @export
cov.taper <- function(cov, supLen=1) {
  function(X, beta, ...) {
    # Formals for 'cov'
    args = list(..., X=X, beta=beta)

    # Sparse distance matrix; don't need full K
    D = distanceMatrix(X, max.dist=supLen)
    args$D = D
    K = do.call(cov, args)

    # Sparse tapered K
    taperMat = cov.MR(X, beta=log(c(1, rep(supLen, ncol(X)))))
    taperMat * K
  }
}

#' Taper a covariance function (partial derivatives)
#'
#' Works as cov.taper, except this function takes a function producing the
#' partial derivatives of a covariance function. This modifies the function
#' to produce the partial derivatives of the tapered covariance function.
#'
#' @param cov.d The partial derivative function of the covariance function to
#'    be tapered
#' @param supLen The support length of the MR
#' @export
cov.taper.d <- function(cov.d, supLen=1) {
  function(X, beta, ...) {
    # Formals for 'cov.d'
    args = list(..., X=X, beta=beta)

    # Sparse distance matrix so we don't compute full dKs
    D = distanceMatrix(X, max.dist=supLen)
    args$D = D
    dK = do.call(cov.d, args)

    # Apply taper to supported elements of each dK
    taperMat = cov.MR(X, beta=log(c(1, rep(supLen, ncol(X)))))
    lapply(dK, function(dKi) dKi * taperMat)
  }
}

cov.triangular <- function(X, beta, ...) {
  stopifnot(is.matrix(X))
  D   = distanceMatrix(X, max.dist=exp(beta[2]))
  D@x = exp(beta[1])*(1 - D@x/exp(beta[2]))
  return(D)
}

cov.triangular.d <- function(X, beta, ...) {
  stopifnot(is.matrix(X))
  D   = distanceMatrix(X, max.dist=exp(beta[2]))

  dKdBeta1 = D
  dKdBeta1@x = 1 - D@x/exp(beta[2])

  dKdBeta2 = D
  dKdBeta2@x = exp(beta[1])*D@x/(exp(beta[2])^2)
  return(list(dKdBeta1, dKdBeta2))
}
