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
#' @import maxLik
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

  ## BIC
  n   = nrow(stpcaObj$Xc) # Number of samples
  d   = ncol(stpcaObj$Xc) # Original dimensionality
  k   = ncol(stpcaObj$W) # Latent dimensionality
  dof = d*k - 0.5*k*(k-1) + 3 # Degrees of Freedom for PPCA
  stpcaObj$bic = -2*stpcaObj$ll + dof*log(n)

  return(stpcaObj)
}

#' Update theta to be the maximum-a-posteriori value using Expectation
#' Maximisation.
#'
#' @param stpcaObj an stpca object.
#' @param maxit.inner number of EM iterations
#' @return An \code{stpca} object.
#' @export
#' @include statistical-quantities.R
#' @include EM.R
stpca.iterate.theta <- function(stpcaObj, maxit.inner=10) {
  vars = within(unclass(stpcaObj), {
    for (iteration in seq_len(maxit.inner)) {
      # Expectation step
      expectations = EM.E(Xc, W, sigSq)
      V = expectations$V; Vvar = expectations$Vvar

      # Maximization step for W
      W = EM.M.W(Xc, sigSq, V, Vvar, K)

      # Second expectation Step
      expectations = EM.E(Xc, W, sigSq)
      V = expectations$V; Vvar = expectations$Vvar

      # Maximization step for sigma^2
      sigSq = EM.M.sigSq(Xc, W, V, Vvar)

      # Remove nonidentifiability VW^T = (VR)(WR)^T by setting R=I
      # This also has the effect of making each column of W an eigenvector of
      # cov[ X | \beta, \sigma^2 ]
      W.svd = svd(W)
      W     = W %*% W.svd$v
      V     = V %*% W.svd$v

      # Calculate new log posterior
      ll = stpca.log_likelihood(Xc, W, rep(0,ncol(Xc)), sigSq)
      lp = stpca.log_posterior(Xc, K, W, rep(0,ncol(Xc)), sigSq)
      lps[length(lps)+1] = lp
    }
  })

  stpcaObj$W     = vars$W
  stpcaObj$sigSq = vars$sigSq
  stpcaObj$V     = vars$V
  stpcaObj$Vvar  = vars$Vvar
  stpcaObj$ll    = vars$ll
  stpcaObj$lp    = vars$lp
  stpcaObj$lps   = vars$lps
  return(stpcaObj)
}

#' Update beta to be the maximum marginal likelihood value using BFGS.
#'
#' @param stpcaObj an stpca object.
#' @return An \code{stpca} object.
#' @export
#' @include statistical-quantities.R
#' @include evidence-approximation.R
#' @import maxLik
stpca.iterate.beta <- function(stpcaObj) {
  vars = within(unclass(stpcaObj), {
    # Set up objective and gradient functions
    logLik = function(beta_) {
      K = covar.fn(locations, beta=beta_)
      stpca.log_evidence(Xc, K, W, rep(0, d), sigSq)
    }

    grad = function(beta_) {
      K  = covar.fn(  locations, beta=beta_)
      dK = covar.fn.d(locations, beta=beta_)
      stpca.log_evidence_d(Xc, K, W, rep(0, d), sigSq, beta_, dK)
    }

    # Do optimisation & get maxlik beta
    optObj = maxBFGS(logLik, grad, start=beta, constraints=constraints)
    beta   = optObj$estimate

    # Add 95% confidence intervals for beta
    stdErr = sqrt(-diag(solve(maxLik::hessian(optObj))))
    ci95 = matrix(NA, nrow=2, ncol=length(beta),
                  dimnames=list(c("upper", "lower"), names(beta)))
    ci95["upper",] = beta + stdErr
    ci95["lower",] = beta - stdErr
    attr(beta, "ci95") = ci95

    # Recompute values depending on beta
    K    = covar.fn(locations, beta=beta, D=D, max.dist=max.dist)
    H    = stpca.H(Xc, W, rep(0,ncol(Xc)), sigSq, K)
    log_evidence = optObj$maximum
    log_evidences = c(log_evidences, log_evidence)
  })

  stpcaObj$H    = vars$H
  stpcaObj$K    = vars$K
  stpcaObj$beta = vars$beta
  stpcaObj$log_evidence  = vars$log_evidence
  stpcaObj$log_evidences = vars$log_evidences
  return(stpcaObj)
}
