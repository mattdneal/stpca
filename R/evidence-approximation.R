#' Compute the laplace approximation to the log evidence given the MAP
#' parameters K, mu, sigSq as well as the prior covariance matrix K.
#' Note that this is multiplied by an UN-KNOWN CONSTANT due to the flat
#' priors over mu and sigSq. However, this unknown constant is always
#' the same regardless of k and K, so this may be used to compute
#' meaningful bayes factors between SPCA models.
#'
#' @param X Data
#' @param K Prior covariance matrix
#' @param W Loadings matrix
#' @param mu
#' @param sigSq
#' @return Approximate log evidence
#' @export
spca.log_evidence <- function(X, K, W, mu, sigSq) {
  if (is(X, "spca")) {
    spcaObj = X
    X = spcaObj$X
    K = spcaObj$K
    W = spcaObj$W
    mu = spcaObj$mu
    sigSq = spcaObj$sigSq
  }

  n = nrow(X)
  d = ncol(X)
  k = ncol(W)

  # Centered X
  Xc = sweep(X, 2, mu)

  # Compute C^{-1}, which is used all over the place
  R = Matrix::chol(crossprod(W) + sigSq*diag(k))
  Cinv = Matrix(Diagonal(d) - Matrix(crossprod(forwardsolve(t(R), t(W)))))/sigSq

  logDetH = 0

  # Compute each of the blocks of H corresponding to each w_i, and the log
  # determinant of this block to the comulative log determinant
  #for (i in 1:k) {
  #  Hw = spca.H.W(X, W[,i], mu, K)
  #  logDetH = logDetH + as.numeric(determinant(Hw, logarithm=TRUE)$modulus)
  #}
  HW = spca.H.W(X, W, mu, sigSq, K)
  logDetH = sum(vapply(HW, function(Hw) {
     as.numeric(determinant(Hw, logarithm=TRUE)$modulus)
  }, numeric(1)))

  # Compute the mu block of H, add the log det to the cumulative total
  HmuLogDet = as.numeric(determinant(n*Cinv, logarithm=TRUE)$modulus)
  logDetH   = logDetH + HmuLogDet

  # Compute the sigSq block of H & add log det to cumulative total
  tmp = (sum(diag(Xc%*%Cinv%*%Cinv%*%Cinv%*%t(Xc))) -
                     0.5*n*sum(diag(Cinv%*%Cinv)))
  if (tmp<=0) {browser()}
  HsigSqLogDet = log(sum(diag(Xc%*%Cinv%*%Cinv%*%Cinv%*%t(Xc))) -
                     0.5*n*sum(diag(Cinv%*%Cinv)))
  logDetH      = logDetH + HsigSqLogDet

  # Laplace-approximated log evidence
  logZ = (spca.log_posterior(X, K, W, mu, sigSq) +
          (0.5*(d*k+d+1))*log(2*pi) -
          0.5*logDetH)
  return(logZ)
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
spca.log_bayes_factor <- function(X, K1, W1, mu1, sigSq1, K2, W2, mu2, sigSq2) {
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

  ev1 = spca.log_evidence(X1, K1, W1, mu1, sigSq1)
  ev2 = spca.log_evidence(X2, K2, W2, mu2, sigSq2)
  return(ev1-ev2)
}

#' Internal function for computing the w_i block of H
#'
#' @param X Data
#' @param W Loadings matrix
#' @param mu
#' @param sigSq
#' @param K Prior covariance matrix
#' @return H_{w_i}
#' @import Matrix
#' @examples
#' d=10; k=3; n=10
#' X = matrix(rnorm(n*d), ncol=d)
#' W = matrix(rnorm(d*k), ncol=k)
#' mu = rnorm(d)
#' sigSq = rnorm(1)^2
#' K = cov.SE(matrix(1:10, ncol=1), beta=log(c(1, 3)))
#'
#' Hw1.analytic = spca.H.W(X, W, mu, sigSq, K)[[1]]
#' Hw1.numeric  = Matrix(hessian(function(w) {
#'   W_ = W
#'   W_[,1] = w
#'   -spca.log_posterior(X, K, W_, mu, sigSq)
#' }, x=W[,1]))
#'
#' stopifnot(all.equal(HW.analytic, HW.numeric))
spca.H.W <- function(X, W, mu, sigSq, K) {
  n = nrow(X)
  d = ncol(X)
  k = ncol(W)
  Xc = sweep(X, 2, mu)

  R = Matrix::chol(crossprod(W) + sigSq*diag(k))
  Cinv = Matrix(Diagonal(d) - Matrix(crossprod(forwardsolve(t(R), t(W)))))/sigSq

  HW = list()
  for (k_ in 1:k) {
    wi = W[,k_,drop=F]
    HW[[k_]] = Matrix(
      solve(K) +
      Cinv*as.numeric(t(wi)%*%Cinv%*%t(Xc)%*%Xc%*%Cinv%*%wi -
                      n*t(wi)%*%Cinv%*%wi + n) +
      Cinv%*%(
        t(Xc)%*%Xc%*%Cinv%*%wi%*%t(wi) +
        wi%*%t(wi)%*%Cinv%*%t(Xc)%*%Xc +
        as.numeric(t(wi)%*%Cinv%*%wi - 1)*t(Xc)%*%Xc -
        n*wi%*%t(wi)
      )%*%Cinv, forceCheck=TRUE) # TODO: use crossprod on wi%*%t(wi)
  }

  return(HW)
}

min.f.generator <- function(X, W, mu, sigSq, locations, covar.fn, covar.fn.d, D=NA, max.dist=Inf, sparse=FALSE) {

  n = nrow(X)
  d = ncol(X)
  k = ncol(W)

  min.f.sparse = function(beta_) {
    stop("Not yet implemented!")
  }

  min.f.dense = function(beta_) {
    K_  = covar.fn(locations, beta=beta_, D=D, max.dist=max.dist)
    dK_ = covar.fn.d(locations, beta=beta_, D=D, max.dist=max.dist)

    if (any(is.na(K_@x)) | any(is.nan(K_@x)) | any(is.infinite(K_@x))) { browser() }

    # TODO: Calculate determinant from decomposition
    logDetK = as.numeric(determinant(K_, logarithm=TRUE)$modulus)
    trWtKinvW = sum(diag(crossprod(W, solve(K_, W))))

    ## Terms used in f and gradient calcs: HwSum and logDetHwSum
    HW = spca.H.W(X, W, mu, sigSq, K_)
    HwSum = Reduce('+', HW)
    logDetHwSum = vapply(HW, function(Hw) {
      as.numeric(determinant(Hw, logarithm=TRUE)$modulus)
    }, numeric(1))
    grTerm3 = solve(K_, t(solve(K_, t(HwSum))))

    ## Function value
    f.val = 0.5*(k*logDetK + trWtKinvW + logDetHwSum)

    ## Function gradients
    f.gr = numeric(length(dK_))
    for (i in 1:length(f.gr)) {
      grTerm1 = k * sum(diag(solve(K_, dK_[[i]])))
      KinvW = solve(K_, W)
      grTerm2 = sum(vapply(1:k, function(k_) {
        as.numeric(KinvW[,k_] %*% dK_[[i]] %*% KinvW[,k_])
      }, numeric(1)))
      f.gr[i] = 0.5*(grTerm1 - grTerm2 - sum(diag(grTerm3%*%dK_[[i]]))) # TODO: This might need a minus sign
    }

    return(list(f=f.val, df=f.gr))
  }


  if(sparse) {
    return(min.f.sparse)
  } else {
    return(min.f.dense)
  }
}
