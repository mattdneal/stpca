#' Compute all the blocks of H.
#'
#' @param X Data
#' @param W Loadings matrix
#' @param mu
#' @param sigSq
#' @param K Prior covariance matrix
#' @return H
#' @import Matrix
#' @examples
#' set.seed(1)
#' d=10; k=3; n=1000
#' X = matrix(rnorm(n*d), ncol=d)
#' W = matrix(rnorm(d*k), ncol=k)
#' mu = rnorm(d)
#' sigSq = rnorm(1)^2
#' K = cov.SE(matrix(1:10, ncol=1), beta=log(c(2, 3)))
#'
#' library(numDeriv)
#' library(Matrix)
#'
#' #Test that the analytic hessian for mu & sigSq matches numerical Hessian.
#' H.analytic = stpca:::compute_H(X, W, mu, sigSq, K)
#' HsigSq.numeric = Matrix(numDeriv::hessian(function(sigSq_) {
#'   -log_posterior(X, K, W, mu, sigSq_)
#' }, x=sigSq))
#' stopifnot(all.equal(H.analytic$sigSq, HsigSq.numeric,
#'                     tolerance=1e-8))
#'
#' Hmu.numeric = Matrix(numDeriv::hessian(function(mu_) {
#'   -log_posterior(X, K, W, mu_, sigSq)
#' }, x=mu))
#' stopifnot(all.equal(H.analytic$mu, Hmu.numeric, tolerance=1e-6))
compute_H <- function(X, W, mu, sigSq, K) {
  n = nrow(X)
  d = ncol(X)
  k = ncol(W)
  Xc = sweep(X, 2, mu)

  R = Matrix::chol(crossprod(W) + sigSq*diag(k))
  Cinv = Matrix(Diagonal(d) - Matrix(crossprod(forwardsolve(t(R), t(W)))))/sigSq

  HW  = compute_H_W(X, W, mu, sigSq, K)

  Hmu = n*Cinv

  # TODO: Can definitely optimise this
  HsigSq = Matrix(sum(diag(Xc%*%Cinv%*%Cinv%*%Cinv%*%t(Xc))) -
                     0.5*n*sum(diag(Cinv%*%Cinv)))
  H = list()
  H[paste("w",1:length(HW),sep='')] = HW
  H["mu"]    = Hmu
  H["sigSq"] = HsigSq
  return(H)
}

#' Compute all the w_i blocks of H
#'
#' @param X Data
#' @param W Loadings matrix
#' @param mu
#' @param sigSq
#' @param K Prior covariance matrix
#' @return H_{w_i}
#' @import Matrix
#' @examples
#' set.seed(1)
#' d=10; k=3; n=10
#' X = matrix(rnorm(n*d), ncol=d)
#' W = matrix(rnorm(d*k), ncol=k)
#' mu = rnorm(d)
#' sigSq = rnorm(1)^2
#' K = cov.SE(matrix(1:10, ncol=1), beta=log(c(1, 3)))
#'
#' library(numDeriv)
#' library(Matrix)
#' Hw1.analytic = stpca:::compute_H_W(X, W, mu, sigSq, K)[[1]]
#' Hw1.numeric  = Matrix(numDeriv::hessian(function(w) {
#'   W_ = W
#'   W_[,1] = w
#'   -log_posterior(X, K, W_, mu, sigSq)
#' }, x=W[,1]))
#'
#' stopifnot(all.equal(Hw1.analytic, Hw1.numeric))
compute_H_W <- function(X, W, mu, sigSq, K) {
  n = nrow(X)
  d = ncol(X)
  k = ncol(W)
  Xc = Matrix(sweep(X, 2, mu))

  M = Matrix(crossprod(W) + sigSq*diag(k))
  R = Matrix::chol(M)
  Cinv = (Diagonal(d) - Matrix(crossprod(forwardsolve(t(R), t(W)))))/sigSq

  MinvWt = solve(M, t(W))

  invSuccess=FALSE
  try({
    Kinv = solve(K, sparse=FALSE)
    invSuccess=TRUE
  }, silent=TRUE)
  if (!invSuccess) { stop("Could not invert K") }
  Kinv = forceSymmetric(Kinv) # Kinv isn't symmetric here; failure of solve method?

  HW = list()
  for (k_ in 1:k) {
    wi = Matrix(W[,k_,drop=F])

    wtCinvw = crossprod(wi, Cinv%*%wi)
    term2 = Cinv*as.numeric(crossprod(Xc%*%(Cinv%*%wi))
                            - n*wtCinvw + n)

    # These 3 (d*d) matrices have all been kept symmetric to reduce computation
    term3A = 2*symmpart(tcrossprod(crossprod(Xc, Xc %*% (Cinv %*% wi)), wi))
    term3B = as.numeric(wtCinvw - 1)*crossprod(Xc)
    term3C = -n*tcrossprod(wi)
    ABCsum = term3A + term3B + term3C

    AW = ABCsum %*% W
    WtAW = crossprod(W, AW)
    AWMinvWt = AW%*%MinvWt
    term3 = (
      ABCsum -
      2*symmpart(AWMinvWt) +
      forceSymmetric(crossprod(MinvWt, WtAW) %*% MinvWt)
    )/(sigSq*sigSq)

    Hwk = Kinv + term2 + term3
    HW[[k_]] = Hwk
  }

  return(HW)
}

#' Compute the partial derivatives of log(det(H)) with respect to the
#' hyperparameters with the value of beta provided.
#'
#' @param K Prior covariance matrix
#' @param KD Prior covariance matrix derivatives
#' @param HW list of blocks H_{w_i}
#' @return Partial derivatives of log|H|
log_det_H_d <- function(K, KD, HW) {
  logDetH.d = numeric(length(KD))
  for (i in seq_along(KD)) {
    logDetH.d[i] = -sum(vapply(HW, function(Hw) {
      sum(diag( solve(K%*%Hw%*%K, KD[[i]]) ))
    }, numeric(1)))
  }
  return(logDetH.d)
}
