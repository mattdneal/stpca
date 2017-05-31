#' Performs SPCA
#'
#' @param X Data
#' @param k Latent dimensionality
#' @param locations the coordinates of each dimenion in X
#' @param covar.fn covariance function to generate K_beta
#' @param covar.fn.d gradient of the covariance function with respect to hyperparameters
#' @param beta0 initial hyperparameters
#' @param trace amount of reporting. 0=none, 1=low, 2=high
#' @param report_iter Number of iterations between reports
#' @param max.dist Maximum distance between features to consider
#' @param maxit number of inner iterations
#' @param maxit.outer number of outer iterations
#' @return An \code{spca} object.
#' @export
#' @include util.R
#' @include evidence-approximation.R
#' @import gsl
#' @import fields
#' @examples
#' library(fields)
#' data(ozone2)
#'
#' # Missing data: Replace missing values by column means
#' X = ozone2$y
#' for (col in 1:ncol(X)) {
#'   ind = is.na(X[,col])
#'   X[ind,col] = mean(X[,col], na.rm=TRUE)
#' }
#' X = X/sd(X) # Scale for numerical reasons
#'
#' locations = ozone2$lon.lat
#' locations = apply(locations, 2, function(col) (col-min(col))/(max(col)-min(col)))
#'
#' model.spca = spca(X, 3, locations, cov.SE, cov.SE.d, beta0=log(c(1, 0.5)),
#'                   maxit=20, maxit.outer=3, trace=0)
spca <- function(X, k, locations, covar.fn, covar.fn.d=NULL, beta0=c(),
                 trace=0, report_iter=10, max.dist=Inf,
                 maxit=20, maxit.outer=5) {

  ## Define commonly used variables.
  Xc = scale(X, scale=FALSE) # Centered data: nxd
  mu = attr(Xc, "scaled:center")
  n  = nrow(X)               # Number of samples
  d  = ncol(X)               # Original dimensionality

  ## Perform sanity checks.
  stopifnot(ncol(X) > k) # Cannot deal with complete/overcomplete case
  stopifnot(nrow(X) >= k) # TODO: Check if I can deal with equality case

  ## Initialize W and sigma^2 from PPCA
  covar.svd = svd(Xc/sqrt(n), nu=0, nv=k)
  covar.eigval = covar.svd$d^2
  sigSq = sum(covar.eigval[-(1:k)])/(d-k)
  W     = covar.svd$v %*% diag(sqrt(covar.eigval[1:k] - sigSq), ncol=k, nrow=k)

  if (sigSq < 1e-10) {warning("The data provided lie close to a subspace of ",
    "dimensionality equal to or lower than the k provided; spca may fail due ",
    "to producing a degenerate probability model. Maybe pick a smaller k?")}

  if (trace>=1) {
    print(paste("Starting spca with", length(beta0), "hyperparameters"))
  }

  beta = beta0
  D    = distanceMatrix(locations, max.dist=max.dist)
  if (length(beta)>0) {
    K    = covar.fn(locations, beta=beta, D=D, max.dist=max.dist)
  } else { # No HPs
    K    = covar.fn(locations, D=D, max.dist=max.dist)
  }
  stopifnot(is(K, "Matrix"))

  lp   = spca.log_posterior(X, K, W, mu, sigSq) # Current log posterior
  lps  = numeric(maxit) # Record of log posteriors for monitoring convergence

  outerConverged = FALSE
  innerConverged = FALSE
  iteration = 0
  outerIteration = 0
  while (!outerConverged) {
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
      if (all(WtW==0)) {stop(paste("SPCA has failed due to numerical instability.",
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
      lpNew = spca.log_posterior(X, K, W, mu, sigSq)
      if (iteration >= maxit) { # Check for maximum iterations reached. If so, print.
        if (trace>0) {
          print(paste("Convergence criteria reached:", iteration, "iterations"))
        }
        innerConverged=TRUE
      } else if (trace==1 && (iteration%%report_iter)==0) {
        print(paste("Iteration ", iteration, ": log likelihood = ",
                    round(lpNew, 4), " (increase=", round(lpNew-lp, 4),")", sep=''))
      }

      lp = lpNew
      iteration = iteration + 1
      lps[iteration] = lp
    } # end 'innerConverged' loop

    ## Remove nonidentifiability VW^T = (VR)(WR)^T by setting R=I
    ## This also has the effect of making each column of W an eigenvector of cov[ X | \beta, \sigma^2 ]
    W.svd = svd(W)
    W     = W.svd$u %*% diag(W.svd$d, nrow=k, ncol=k)
    E_V1  = E_V1 %*% W.svd$v

    # TODO: Move out to a different function & get unit tests. Same w/ 'inner' EM loop.
    ########################
    ## Tune Hyperparameters
    ########################
    outerConverged = (outerIteration>=maxit.outer)
    if (!outerConverged & (length(beta0) > 0)) { # There are HPs to tune
      evidence = spca.log_evidence(X, K, W, mu, sigSq)
      if (trace>=1) {
        print(paste("Outer iteration ", outerIteration, ": log evidence=",
                    round(evidence, 4), sep=''))
      }

      min.f = function(beta_) {
        K_ = covar.fn(locations, beta=beta_, D=D, max.dist=max.dist)
        if (any(is.na(K_@x)) | any(is.nan(K_@x)) | any(is.infinite(K_@x))) { browser() }
        return(-spca.log_evidence(X, K_, W, mu, sigSq))
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

  if (iteration < maxit) { # Trim unused lps entries
    lps = lps[1:iteration]
  }

  # TODO: Make largest abs value positive
  # Identify directionality of each component by fixing sign of 1st element to be +ve
  #P = diag(sign(W[1,]), nrow=k, ncol=k)
  #W = W %*% P
  #E_V1 = E_V1 %*% P

  ll = spca.log_likelihood(X, W, mu, sigSq)
  dof = d*k - 0.5*k*(k-1) + 3 + length(beta0) # Degrees of Freedom for PPCA + #HPs
  bic = -2*ll + dof*log(n)

  H = spca.H(X, W, mu, sigSq, K)

  spcaObj = list(X     = X,
                 W     = W,
                 sigSq = sigSq,
                 mu    = mu,
                 V     = E_V1,
                 ll    = ll,
                 lp    = lp,
                 lps   = lps,
                 bic   = bic,
                 beta  = beta,
                 K     = K,
                 H     = H,
                 covar.fn = covar.fn,
                 locations = locations,
                 covar.fn.d = covar.fn.d,
                 log_evidence = evidence)
  class(spcaObj) = "spca"
  return(spcaObj)
}

#' Does not yet work with 'new' SPCA architecture
spca.continue <- function(spcaObj, trace=0, report_iter=10, max.dist=Inf,
                          maxit=10, tol=1e-2, ucminf.control=list()) {
  stop("Does not yet work with 'new' SPCA architecture")
  newPrcaObj = spca(X=spcaObj$X, k=ncol(spcaObj$W), locations=spcaObj$locations,
                    covar.fn=spcaObj$covar.fn, covar.fn.d=spcaObj$covar.fn.d,
                    beta0=spcaObj$beta, trace=trace, report_iter=report_iter,
                    max.dist=max.dist, maxit=maxit, tol=tol,
                    ucminf.control=ucminf.control)
  return(newPrcaObj)
}

#' Calculate the log likelihood for SPCA with given parameters
#'
#' @param X Data
#' @param W
#' @param mu
#' @param sigSq
#' @return log likelihood (numeric)
#' @export
spca.log_likelihood <- function(X, W, mu, sigSq) {
  if (is.vector(X)) {X = matrix(X, nrow=1)}

  Xc = sweep(X, 2, mu)

  d = nrow(W)
  k = ncol(W)
  n = nrow(X)

  W.sv = svd(W, nu=0, nv=0)$d
  R    = chol(crossprod(W) + sigSq*diag(k))

  const     = n*d*log(2*pi)
  nLogDetC  = n*((d-k)*log(sigSq) + sum(log(W.sv^2 + sigSq)))
  trXCinvXt = (norm(Xc, 'F')^2 -
               norm(forwardsolve(t(R), t(W))%*%t(Xc), 'F')^2)/sigSq
  return(-0.5*(const + nLogDetC + trXCinvXt))
}

#' Calculate the *un-normalised* log prior for SPCA with given loadings matrix
#'
#' @param K Prior covariance matrix
#' @param W Loadings matrix
#' @return un-normalised log prior (numeric)
#' @export
spca.log_prior <- function(K, W) {
  d = nrow(W)
  k = ncol(W)

  # This special case for sparse matrices is more numerically stable
  if (is(K, "sparseMatrix")) {
    if (all(K@x==0) | all(K@x==Inf)) {
      return(-Inf)
    }

    Kc = Cholesky(K, LDL=TRUE, pivot=TRUE)

    # Fast calculation of the log determinant of K
    #logDetK = as.numeric(-determinant(solve(Kc, system='D'), log=TRUE)$modulus)

    logDetK = as.numeric(determinant(K, logarithm=TRUE)$modulus)
    if (logDetK==Inf) {return(-Inf)}

    KinvW     = solve(Kc, W, system="A")
    trWtKinvW = sum(vapply(1:k, function(k_) (W[,k_] %*% KinvW[,k_])[1,1],
                           numeric(1)))
  } else {
    K = Matrix(K)
    # This is more stable than a base matrix solution, but not for sparse
    # Matrices (dealt with above)
    if (all(K@x==0) | all(K@x==Inf)) {return(-Inf)}
    trWtKinvW = sum(diag(Matrix::crossprod(W, solve(K, W))))
    # TODO: Calculate determinant from decomposition
    logDetK = as.numeric(determinant(K, logarithm=TRUE)$modulus)
  }

  return(-0.5*( k*d*log(2*pi) +
                k*logDetK +
                trWtKinvW))
}

#' Calculate the *un-normalised* log posterior for SPCA with given parameters
#'
#' @param X Data
#' @param K Prior covariance matrix
#' @param W Loadings matrix
#' @param mu
#' @param sigSq
#' @return un-normalised log posterior (numeric)
#' @export
spca.log_posterior <- function(X, K, W, mu, sigSq) {
  return(spca.log_likelihood(X, W, mu, sigSq) + spca.log_prior(K, W))
}

#' De-noise a sample using a trained \code{spca} object.
#'
#' @param object An \code{spca} object returned from a call to \code{spca}
#' @param samples samples to de-noise
#' @return De-noised samples of the same dimensionality as the parameter \code{samples}
#' @export
predict.spca <- function(object, samples) {
  if (missing(samples)) {
    return(object$V)
  }

  stopifnot(is.matrix(samples))

  k     = ncol(object$W) # Latent dimensionality
  muMat = matrix(rep(object$mu, nrow(samples)), nrow=nrow(samples), byrow=TRUE)

  # Latent representation of new samples
  Mchol  = chol(crossprod(object$W) + object$sigSq*diag(k))
  latent = t(backsolve(Mchol, forwardsolve(t(Mchol), t((samples-muMat)%*%object$W))))

  # 'samples' projected onto PrCA model
  proj = tcrossprod(latent, object$W) + muMat
  return(proj)
}

spca.simulate <- function(object, n=1) {
  stop("Implement me!")
}
