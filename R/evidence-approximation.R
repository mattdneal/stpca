#' Compute the laplace approximation to the log evidence given the MAP
#' parameters K, mu, sigSq as well as the prior covariance matrix K.
#' Note that this is multiplied by an UN-KNOWN CONSTANT due to the flat
#' priors over mu and sigSq. However, this unknown constant is always
#' the same regardless of k and K, so this may be used to compute
#' meaningful bayes factors between StPCA models.
#'
#' @param X Data
#' @param K Prior covariance matrix
#' @param W Loadings matrix
#' @param mu
#' @param sigSq
#' @return Approximate log evidence
#' @export
stpca.log_evidence <- function(X, K, W, mu, sigSq) {
  if (is(X, "stpca")) {
    stpcaObj = X
    X = stpcaObj$X
    K = stpcaObj$K
    W = stpcaObj$W
    mu = stpcaObj$mu
    sigSq = stpcaObj$sigSq
  }

  if (!any(is.finite(diag(K)))) {
    return(-Inf)
  }

  n = nrow(X)
  d = ncol(X)
  k = ncol(W)

  # If the inversion cannot be done, logZ defaults to -Inf
  logZ = -Inf
  try({
    H = stpca.H(X, W, mu, sigSq, K)
    logDetH = sum(vapply(H, function(Hblock) {
       as.numeric(determinant(Hblock, logarithm=TRUE)$modulus)
    }, numeric(1)))

    # Laplace-approximated log evidence
    logZ = (stpca.log_posterior(X, K, W, mu, sigSq) +
            (0.5*(d*k+d+1))*log(2*pi) -
            0.5*logDetH)
  }, silent=TRUE)
  return(logZ)
}

#' Compute the derivative of the approximate log evidence with respect to the
#' value of beta provided.
#'
#' @param X Data
#' @param K Prior covariance matrix
#' @param W Loadings matrix
#' @param mu
#' @param sigSq
#' @param beta
#' @param dK
#' @return Partial derivatives of approximate log evidence
#' @export
stpca.log_evidence_d <- function(X, K, W, mu, sigSq, beta, dK) {
  success=FALSE
  try({HW = stpca.H.W(X, W, mu, sigSq, K); success=TRUE})

  if(!success) {return(lapply(seq_along(beta), function(b) 0))}

  logPrior.d = stpca.log_prior_d(W, beta, K, dK)
  logDetH.d  = stpca.log_det_H_d(K, dK, HW)

  logEvidence.d = vapply(seq_along(beta), function(i) {
    logPrior.d[i] - 0.5*logDetH.d[i]
  }, numeric(1))
  names(logEvidence.d) = names(beta)
  return(logEvidence.d)
}

#' Compute the partial derivatives of the log prior with respect to the
#' hyperparameters with the value of beta provided.
#'
#' @param W Loadings matrix
#' @param beta hyperparameter values
#' @param K Prior covariance matrix
#' @param dK Prior covariance matrix derivatives
#' @return Partial derivatives of log prior
#' @export
stpca.log_prior_d <- function(W, beta, K, dK) {
  deriv = numeric(length(beta))
  KinvW = solve(K, W)

  for (i in seq_along(beta)) {
    term1 = ncol(W) * sum(diag( solve(K, dK[[i]])  ))

    term2 = -sum(diag( tcrossprod(KinvW)%*%dK[[i]] ))
    #term2b = -sum(vapply(1:ncol(W), function(k_) {
    #  as.numeric(KinvW[,k_] %*% dK[[i]] %*% KinvW[,k_])
    #}, numeric(1)))

    deriv[i] = -0.5*(term1 + term2)
  }
  names(deriv) = names(beta)
  return(deriv)
}

#' Compute the partial derivatives of log(det(H)) with respect to the
#' hyperparameters with the value of beta provided.
#'
#' @param K Prior covariance matrix
#' @param dK Prior covariance matrix derivatives
#' @param HW list of blocks H_{w_i}
#' @return Partial derivatives of log|H|
#' @export
stpca.log_det_H_d <- function(K, dK, HW) {
  logDetH.d = numeric(length(dK))
  for (i in seq_along(dK)) {
    logDetH.d[i] = -sum(vapply(HW, function(Hw) {
      sum(diag( solve(K%*%Hw%*%K, dK[[i]]) ))
    }, numeric(1)))
  }
  return(logDetH.d)
}

#' Compute bayes factor
#'
#' @param X Data
#' @param K1
#' @param W1
#' @param mu1
#' @param sigSq1
#' @param K2
#' @param W2
#' @param mu2
#' @param sigSq2
#' @return Log bayes factor; model 1 is numerator
#' @export
stpca.log_bayes_factor <- function(X, K1, W1, mu1, sigSq1, K2, W2, mu2, sigSq2) {
  only2argsspecified = (!missing(X) & !missing(K1) & missing(W1) & missing(mu1)
                        & missing(sigSq1) & missing(K2) & missing(W2)
                        & missing(mu2) & missing(sigSq2))

  if(only2argsspecified) {
    model1 = X
    model2 = K1

    X1     = model1$X
    K1     = model1$K
    W1     = model1$W
    mu1    = model1$mu
    sigSq1 = model1$sigSq

    X2     = model2$X
    K2     = model2$K
    W2     = model2$W
    mu2    = model2$mu
    sigSq2 = model2$sigSq
  }

  ev1 = stpca.log_evidence(X1, K1, W1, mu1, sigSq1)
  ev2 = stpca.log_evidence(X2, K2, W2, mu2, sigSq2)
  return(ev1-ev2)
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
#' Hw1.analytic = stpca:::stpca.H.W(X, W, mu, sigSq, K)[[1]]
#' Hw1.numeric  = Matrix(numDeriv::hessian(function(w) {
#'   W_ = W
#'   W_[,1] = w
#'   -stpca.log_posterior(X, K, W_, mu, sigSq)
#' }, x=W[,1]))
#'
#' stopifnot(all.equal(Hw1.analytic, Hw1.numeric))
stpca.H.W <- function(X, W, mu, sigSq, K) {
  n = nrow(X)
  d = ncol(X)
  k = ncol(W)
  Xc = sweep(X, 2, mu)

  R = Matrix::chol(crossprod(W) + sigSq*diag(k))
  Cinv = forceSymmetric(Diagonal(d) - Matrix(crossprod(forwardsolve(t(R), t(W)))))/sigSq

  HW = list()
  for (k_ in 1:k) {
    wi = W[,k_,drop=F]

    invSuccess=FALSE
    try({
      Kinv = solve(K)
      invSuccess=TRUE
    }, silent=TRUE)
    if (!invSuccess) { stop("Could not invert K") }

    Hwk = Kinv +
      Cinv*as.numeric(tcrossprod(t(wi)%*%Cinv%*%t(Xc))
                      - n*t(wi)%*%Cinv%*%wi + n) +
      tcrossprod(Cinv, crossprod(Cinv, (
        crossprod(Xc)%*%Cinv%*%tcrossprod(wi) +
        tcrossprod(wi)%*%Cinv%*%crossprod(Xc) +
        as.numeric(t(wi)%*%Cinv%*%wi - 1)*crossprod(Xc) -
        n*tcrossprod(wi))))

    HW[[k_]] = 0.5*forceSymmetric(Hwk + t(Hwk))
  }
  return(HW)
}

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
#' H.analytic = stpca:::stpca.H(X, W, mu, sigSq, K)
#' HsigSq.numeric = Matrix(numDeriv::hessian(function(sigSq_) {
#'   -stpca.log_posterior(X, K, W, mu, sigSq_)
#' }, x=sigSq))
#' stopifnot(all.equal(H.analytic$sigSq, HsigSq.numeric,
#'                     tolerance=1e-8))
#'
#' Hmu.numeric = Matrix(numDeriv::hessian(function(mu_) {
#'   -stpca.log_posterior(X, K, W, mu_, sigSq)
#' }, x=mu))
#' stopifnot(all.equal(H.analytic$mu, Hmu.numeric, tolerance=1e-6))
stpca.H <- function(X, W, mu, sigSq, K) {
  n = nrow(X)
  d = ncol(X)
  k = ncol(W)
  Xc = sweep(X, 2, mu)

  R = Matrix::chol(crossprod(W) + sigSq*diag(k))
  Cinv = Matrix(Diagonal(d) - Matrix(crossprod(forwardsolve(t(R), t(W)))))/sigSq

  HW  = stpca.H.W(X, W, mu, sigSq, K)

  Hmu = n*Cinv
  HsigSq = Matrix(sum(diag(Xc%*%Cinv%*%Cinv%*%Cinv%*%t(Xc))) -
                     0.5*n*sum(diag(Cinv%*%Cinv)))
  H = list()
  H[paste("w",1:length(HW),sep='')] = HW
  H["mu"]    = Hmu
  H["sigSq"] = HsigSq
  return(H)
}
