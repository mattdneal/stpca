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

  n = nrow(X)
  d = ncol(X)
  k = ncol(W)

  # Centered X
  Xc = sweep(X, 2, mu)

  # Compute C^{-1}, which is used all over the place
  R = Matrix::chol(crossprod(W) + sigSq*diag(k))
  Cinv = Matrix(Diagonal(d) - Matrix(crossprod(forwardsolve(t(R), t(W)))))/sigSq

  H = stpca.H(X, W, mu, sigSq, K)
  logDetH = sum(vapply(H, function(Hblock) {
     as.numeric(determinant(Hblock, logarithm=TRUE)$modulus)
  }, numeric(1)))

  # Laplace-approximated log evidence
  logZ = (stpca.log_posterior(X, K, W, mu, sigSq) +
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
#' Hw1.numeric  = Matrix(hessian(function(w) {
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
#' HsigSq.numeric = Matrix(hessian(function(sigSq_) {
#'   -stpca.log_posterior(X, K, W, mu, sigSq_)
#' }, x=sigSq))
#' stopifnot(all.equal(H.analytic$sigSq, HsigSq.numeric,
#'                     tolerance=1e-8))
#'
#' Hmu.numeric = Matrix(hessian(function(mu_) {
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

#' Value and gradient of function to be minimized in tuning beta
#'
#' @param X
#' @param W
#' @param mu
#' @param sigSq
#' @param locations
#' @param covar.fn
#' @param covar.fn.d
#' @return A function which, when given a value of beta, returns a value and
#'         gradient of a function to be minimized in beta-tuning.
#' @import Matrix
#' @import numDeriv
#' @examples
#' set.seed(1)
#' n = 100
#' d = 30
#' k = 3
#' X = matrix(rnorm(n*d), ncol=d)
#' W = matrix(rnorm(d*k), nrow=d, ncol=k)
#' mu = rnorm(d)
#' sigSq = rnorm(1)^2
#' locations = matrix(rnorm(d*2), ncol=2)
#' beta = rnorm(2)
#'
#' fdf = stpca:::min.f.generator(X, W, mu, sigSq, locations, cov.SE, cov.SE.d)
#'
#' library(numDeriv)
#' grad.analytic = fdf(beta)$df
#' grad.numeric = grad(function(beta_) {
#'   fdf(beta_)$f
#' }, x=beta)
#' all.equal(grad.analytic, grad.numeric, tolerance=100*.Machine$double.eps^0.5)
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

    HW = stpca.H.W(X, W, mu, sigSq, K_)
    logDetHwSum = sum(vapply(HW, function(Hw) {
      as.numeric(determinant(Hw, logarithm=TRUE)$modulus)
    }, numeric(1)))

    logDetK = as.numeric(determinant(K_, logarithm=TRUE)$modulus)

    TrWtKinvW = sum(diag(Matrix::crossprod(W, solve(K_, W))))

    ## Function value
    f.val = 0.5*(logDetHwSum + k*logDetK + TrWtKinvW)

    ## Function gradients
    f.gr = numeric(length(dK_))
    for (i in 1:length(f.gr)) {
      grTerm1 = -sum(vapply(HW, function(Hw) {
        sum(diag(solve(K_%*%Hw%*%K_, dK_[[i]])))
      }, numeric(1)))

      grTerm2 = k * sum(diag(solve(K_, dK_[[i]])))

      KinvW = solve(K_, W)
      grTerm3 = -sum(vapply(1:k, function(k_) {
        as.numeric(KinvW[,k_] %*% dK_[[i]] %*% KinvW[,k_])
      }, numeric(1)))

      f.gr[i] = 0.5*(grTerm1 + grTerm2 + grTerm3)
    }
    return(list(f=f.val, df=f.gr))
  }

  if(sparse) {
    return(min.f.sparse)
  } else {
    return(min.f.dense)
  }
}
