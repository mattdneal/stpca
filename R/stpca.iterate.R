#' Continue inference in an stpca object.
#'
#' @param stpcaObj an stpca object.
#' @param trace amount of reporting. 0=none, 1=low, 2=high
#' @param report.iter Number of iterations between reports
#' @param max.dist Maximum distance between features to consider
#' @param maxit.inner number of inner iterations
#' @param maxit.outer number of outer iterations
#' @return An \code{stpca} object.
#' @export
#' @include statistical-quantities.R
#' @include evidence-approximation.R
#' @import gsl
#' @import fields
#' @import Matrix
stpca.iterate <- function(stpcaObj, trace=0, report.iter=10,
                          maxit.inner=20, maxit.outer=5) {

  if (length(beta)==0) {
    maxit.outer=0
    if (trace >= 1) {
      print("No hyperparameters; doing no beta-optimisation.")
    }
  }

  for (iteration in seq_len(maxit.outer)) {
    stpcaObj = stpca.iterate.theta(stpcaObj, maxit.inner) # Update theta
    stpcaObj = stpca.iterate.beta(stpcaObj) # Update beta
    if (trace>=1) {
      print(paste("Outer iteration ", iteration, ":\n",
                  "  *log evidence=", round(stpcaObj$evidence, 4),
                  "  *new beta    =", paste(stpcaObj$beta, collapse=','),
                  sep=''))
    }
  }
  stpcaObj = stpca.iterate.theta(stpcaObj, maxit.inner) # Update theta

  ## BIC
  n  = nrow(stpcaObj$X) # Number of samples
  d  = ncol(stpcaObj$X) # Original dimensionality
  k  = ncol(stpcaObj$W) # Latent dimensionality
  dof = d*k - 0.5*k*(k-1) + 3 # Degrees of Freedom for PPCA
  stpcaObj$bic = -2*stpcaObj$ll + dof*log(n)

  return(stpcaObj)
}

stpca.iterate.theta <- function(stpcaObj, maxit.inner=10) {
  stopifnot(all(c("X", "W", "mu", "sigSq", "locations") %in% names(stpcaObj)))
  vars = within(unclass(stpcaObj), {
    for (iteration in seq_len(maxit.inner)) {
      ## Expectation Step
      expectations = expectation(X, W, sigSq)
      E_V1 = expectations$E_V1; E_V2 = expectations$E_V2

      ## Maximization step for sigma^2
      E_V2sum = Reduce('+', E_V2)

      sigSq = (
        norm(X, 'F')^2 -
        2*sum(vapply(1:n, function(n_) E_V1[n_,] %*% t(W) %*% X[n_,], numeric(1))) +
        sum(vapply(1:d, function(d_) W[d_,] %*% E_V2sum %*% W[d_,], numeric(1)))
      )/(n*d)
      sigSq = max(0, sigSq)

      ## Expectation Step
      expectations = expectation(X, W, sigSq)
      E_V1 = expectations$E_V1; E_V2 = expectations$E_V2

      ## Maximization step for W
      xvsum = Reduce('+', lapply(1:n, function(i_) tcrossprod(X[i_,], E_V1[i_,])))
      vvsum.eig = eigen(Reduce('+', E_V2), symmetric=TRUE)
      vvsuminv.eig = list(values=rev(1/vvsum.eig$values),
                          vectors=vvsum.eig$vectors[,k:1])

      C.tilde  = (K %*% xvsum %*% vvsuminv.eig$vectors %*%
                  diag(vvsuminv.eig$values, ncol=k, nrow=k))

      ## I have been unable to speed this up by doing Kc <- Cholesky(K) and
      ## calculating each W.tilde_i with a diagonal update to Kc. Using updown
      ## gives the correct answer but is MUCH slower and using update gives the
      ## wrong answer.
      W.tilde = vapply(1:k, function(i_) {
        if (is(K, "sparseMatrix")) {
          Kc = Cholesky(K, Imult=sigSq*vvsuminv.eig$values[i_], LDL=TRUE, perm=TRUE)
          return(as.vector(Matrix::solve(Kc, C.tilde[,i_], system="A")))
        } else {
          KplusDiag       = K
          Matrix::diag(KplusDiag) = Matrix::diag(KplusDiag) + sigSq*vvsuminv.eig$values[i_]
          return(as.vector(Matrix::solve(KplusDiag, C.tilde[,i_])))
        }
      }, numeric(d))
      W = W.tilde %*% t(vvsuminv.eig$vectors)

      lp = stpca.log_posterior(X, K, W, mu, sigSq)
      lps[length(lps)+1] = lp
    }

    # Remove nonidentifiability VW^T = (VR)(WR)^T by setting R=I
    # This also has the effect of making each column of W an eigenvector of
    # cov[ X | \beta, \sigma^2 ]
    W.svd = svd(W)
    W     = W    %*% W.svd$v
    E_V1  = E_V1 %*% W.svd$v
  })

  stpcaObj$W     = vars$W
  stpcaObj$sigSq = vars$sigSq
  stpcaObj$V     = as.matrix(vars$E_V1)
  stpcaObj$ll    = vars$ll
  stpcaObj$lp    = vars$lp
  stpcaObj$lps   = vars$lps
  return(stpcaObj)
}

stpca.iterate.beta <- function(stpcaObj) {
  stopifnot(all(c("X", "W", "mu", "sigSq", "locations", "covar.fn",
                  "covar.fn.d", "D", "max.dist", "beta") %in% names(stpcaObj)))
  vars = within(unclass(stpcaObj), {
    min.f = function(beta_) {
      K_ = covar.fn(locations, beta=beta_, D=D, max.dist=max.dist)
      if (any(is.na(K_@x)) | any(is.nan(K_@x)) | any(is.infinite(K_@x))) {
        stop(paste("StPCA has failed because one or more elements of K has",
                   "become NA/NaN/Inf"))
      }
      return(-stpca.log_evidence(X, K_, W, mu, sigSq))
    }

    # TODO: Get simultanious evaluationw working; should DOUBLE speed!
    min.fdf = min.f.generator(X, W, mu, sigSq, locations,
                              covar.fn, covar.fn.d, D=D)
    min.f   = function(beta_) {min.fdf(beta_)$f}
    min.df  = function(beta_) {min.fdf(beta_)$df}

    optObj.init = multimin.init(x=beta, #fdf=min.fdf,
                                f=min.f, df=min.df, method="bfgs")
    optObj      = multimin.iterate(optObj.init)

    beta = optObj$x
    K    = covar.fn(locations, beta=beta, D=D, max.dist=max.dist)

    H    = stpca.H(X, W, mu, sigSq, K)

    levidence = stpca.log_evidence(X, K, W, mu, sigSq) # TODO: Remove if this works
    levidence2 = -optObj$f
    print(paste(levidence, "-", levidence2, "=", levidence-levidence2))

    levidence3 = -min.f(beta)
    print(paste(levidence2, "-", levidence3, "=", levidence2-levidence3))
  })

  stpcaObj$H    = vars$H
  stpcaObj$K    = vars$K
  stpcaObj$beta = vars$beta
  stpcaObj$log_evidence = vars$levidence
  return(stpcaObj)
}

expectation <- function(X, W, sigSq) {
  if (all(crossprod(W)==0)) {stop(paste("StPCA has failed due to numerical",
    "instability. Try dividing X by it's largest singular value."))}
  M = Matrix(crossprod(W) + sigSq*diag(ncol(W)))
  E_V1 = t(solve(M, t(X %*% W)))
  E_V2 = lapply(1:nrow(X), function(i_) sigSq*solve(M) + tcrossprod(E_V1[i_,]))
  return(list(E_V1=E_V1, E_V2=E_V2))
}
