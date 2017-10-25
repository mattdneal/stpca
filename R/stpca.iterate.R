#' Continue inference in an stpca object.
#'
#' @param stpcaObj an stpca object.
#' @param trace amount of reporting. 0=none, 1=low, 2=high
#' @param report_iter Number of iterations between reports
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
stpca.iterate <- function(stpcaObj, trace=0, report_iter=10, max.dist=Inf,
                          maxit.inner=20, maxit.outer=5) {

  ## Unpack stpcaObj
  X     = stpcaObj$X
  W     = stpcaObj$W
  sigSq = stpcaObj$sigSq
  mu    = stpcaObj$mu
  beta  = stpcaObj$beta
  D     = stpcaObj$D
  K     = stpcaObj$K
  covar.fn   = stpcaObj$covar.fn
  covar.fn.d = stpcaObj$covar.fn.d
  locations  = stpcaObj$locations

  ## Define commonly used variables.
  Xc = sweep(X, 2, mu) # Centered data
  n  = nrow(X)      # Number of samples
  d  = ncol(X)      # Original dimensionality
  k  = ncol(W)

  lps = c() # Record of log posteriors for monitoring convergence

  outerConverged = FALSE
  outerIteration = 0
  while (!outerConverged) {
    innerConverged = FALSE
    innerIteration = 0
    while (!innerConverged) {
      ##################
      ## EM for sigma^2
      ##################

      ## Expectation Step
      WtW = crossprod(W)
      Minv = chol2inv(chol(WtW + sigSq*diag(k)))
      #Minv2 = Matrix::solve(WtW + sigSq*diag(k)) # TODO: Replace; more stable? # EVEN BETTER: store M; Matrix::solve each system!
      E_V1 = Xc %*% W %*% Minv
      E_V2 = lapply(1:n, function(i_) sigSq*Minv + tcrossprod(E_V1[i_,]))

      ## Maximization step for sigma^2
      E_V2sum = Reduce('+', E_V2)

      sigSq = (
        norm(Xc, 'F')^2 -
        2*sum(vapply(1:n, function(n_) E_V1[n_,] %*% t(W) %*% Xc[n_,], numeric(1))) +
        sum(vapply(1:d, function(d_) W[d_,] %*% E_V2sum %*% W[d_,], numeric(1)))
      )/(n*d)
      sigSq = max(0, sigSq)

      ##################
      ## EM for W
      ##################

      ## Expectation Step
      WtW = crossprod(W)
      if (all(WtW==0)) {stop(paste("StPCA has failed due to numerical instability.",
        "Try dividing X by it's largest singular value to improve numerical",
        "stability"))}
      Minv = chol2inv(chol(WtW + sigSq*diag(k)))
      E_V1 = Xc %*% W %*% Minv
      E_V2 = lapply(1:n, function(i_) sigSq*Minv + tcrossprod(E_V1[i_,]))

      ## Maximization step for W
      xvsum = Reduce('+', lapply(1:n, function(i_) tcrossprod(Xc[i_,], E_V1[i_,])))
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

      ####################################
      ## Convergence criteria & printouts if trace > 0
      ####################################
      innerIteration = innerIteration + 1
      lpNew = stpca.log_posterior(X, K, W, mu, sigSq)
      if (innerIteration >= maxit.inner) { # Check for maximum iterations reached. If so, print.
        if (trace>0) {
          print(paste("Convergence criteria reached:", innerIteration, "iterations"))
        }
        innerConverged=TRUE
      } else if (trace==1 && (innerIteration%%report_iter)==0) {
        print(paste("Iteration ", innerIteration, ": log likelihood = ",
                    round(lpNew, 4), " (increase=", round(lpNew-lp, 4),")", sep=''))
      }

      lp = lpNew
      lps[length(lps)+1] = lp
    } # end 'innerConverged' loop

    # Remove nonidentifiability VW^T = (VR)(WR)^T by setting R=I
    # This also has the effect of making each column of W an eigenvector of
    # cov[ X | \beta, \sigma^2 ]
    W.svd = svd(W)
    W     = W    %*% W.svd$v
    E_V1  = E_V1 %*% W.svd$v

    # TODO: Move out to a different function & get unit tests. Same w/ 'inner' EM loop.
    ########################
    ## Tune Hyperparameters
    ########################
    outerConverged = (outerIteration>=maxit.outer)
    if (!outerConverged & (length(beta) > 0)) { # There are HPs to tune
      levidence = stpca.log_evidence(X, K, W, mu, sigSq)
      if (trace>=1) {
        print(paste("Outer iteration ", outerIteration, ": log evidence=",
                    round(levidence, 4), sep=''))
      }

      min.f = function(beta_) {
        K_ = covar.fn(locations, beta=beta_, D=D, max.dist=max.dist)
        if (any(is.na(K_@x)) | any(is.nan(K_@x)) | any(is.infinite(K_@x))) { browser() }
        return(-stpca.log_evidence(X, K_, W, mu, sigSq))
      }

      #optObj = optimx(par=beta, fn=min.f, method="Nelder-Mead", control=list(
      #  kkt=FALSE, starttests=FALSE, usenumDeriv=TRUE, all.methods=FALSE,
      #  maximize=FALSE, trace=0, dowarn=FALSE
      #))

      ## New version
      min.fdf = min.f.generator(X, W, mu, sigSq, locations, covar.fn, covar.fn.d, D=D)
      min.f   = function(beta_) {min.fdf(beta_)$f}
      min.df  = function(beta_) {min.fdf(beta_)$df}

      #optObj2.init = suppressWarnings(multimin.init(x=beta, fdf=min.fdf,
      #                                f=min.f, df=min.df, method="bfgs"))
      optObj2.init = multimin.init(x=beta, #fdf=min.fdf,
                                   f=min.f, df=min.df, method="bfgs")
      optObj2      = multimin.iterate(optObj2.init)

      #beta = coef(optObj)
      beta = optObj2$x
      K    = covar.fn(locations, beta=beta, D=D, max.dist=max.dist)

      if (trace>=1) { print(paste("  New beta:", paste(beta, collapse=','))) }
      outerIteration = outerIteration+1
    }
  } # end 'outerConverged' loop

  # TODO: Make largest abs value positive
  # Identify directionality of each component by fixing sign of 1st element to be +ve
  #P = diag(sign(W[1,]), nrow=k, ncol=k)
  #W = W %*% P
  #E_V1 = E_V1 %*% P

  ll = stpca.log_likelihood(X, W, mu, sigSq)
  dof = d*k - 0.5*k*(k-1) + 3 + length(beta) # Degrees of Freedom for PPCA + #HPs
  bic = -2*ll + dof*log(n)

  H = stpca.H(X, W, mu, sigSq, K)

  stpcaObj$W     = W
  stpcaObj$sigSq = sigSq
  stpcaObj$V     = as.matrix(E_V1)
  stpcaObj$ll    = ll
  stpcaObj$lp    = lp
  stpcaObj$lps   = c(stpcaObj$lps, lps)
  stpcaObj$bic   = bic
  stpcaObj$beta  = beta
  stpcaObj$D     = D
  stpcaObj$K     = K
  stpcaObj$H     = H
  stpcaObj$log_evidence = levidence
  return(stpcaObj)
}
