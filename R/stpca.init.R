#' Initialises an stpca object.
#'
#' @param X Data
#' @param k Latent dimensionality
#' @param locations the coordinates of each dimenion in X
#' @param covar.fn covariance function to generate K_beta
#' @param covar.fn.d gradient of the covariance function with respect to hyperparameters
#' @param beta0 initial hyperparameters
#' @param trace amount of reporting. 0=none, 1=low, 2=high
#' @param max.dist Maximum distance between features to consider
#' @param constraints constraints on beta, of the form used in maxLik
#' @return An \code{stpca} object.
#' @export
#' @include util.R
#' @include statistical-quantities.R
#' @import Matrix
stpca.init <- function(X, k, locations, covar.fn, covar.fn.d=NULL, beta0=c(),
                       trace=0, max.dist=Inf, constraints=NULL) {

  stopifnot(all(is.finite(X)))

  ## Define commonly used variables.
  Xc  = scale(X, scale=FALSE) # Centered data: nxd
  mu = attr(Xc, "scaled:center")
  n  = nrow(Xc)               # Number of samples
  d  = ncol(Xc)               # Original dimensionality

  ## Perform sanity checks.
  stopifnot(ncol(Xc) > k) # Cannot deal with complete/overcomplete case
  stopifnot(nrow(Xc) >= k) # TODO: Check if I can deal with equality case

  ## Initialize W and sigma^2 from PPCA
  covar.svd = svd(Xc/sqrt(n), nu=0, nv=k)
  covar.eigval = covar.svd$d^2
  sigSq = sum(covar.eigval[-(1:k)])/(d-k)
  W     = covar.svd$v %*% diag(sqrt(covar.eigval[1:k] - sigSq), ncol=k, nrow=k)

  if (sigSq < 1e-10) {warning("The data provided lie close to a subspace of ",
    "dimensionality equal to or lower than the k provided; stpca may fail due ",
    "to producing a degenerate probability model. Maybe pick a smaller k?")}

  if (trace>=1) {
    print(paste("Starting stpca with", length(beta0), "hyperparameters"))
  }

  beta = beta0
  D    = distanceMatrix(locations, max.dist=max.dist)
  if (length(beta)>0) {
    K    = covar.fn(locations, beta=beta, D=D, max.dist=max.dist)
  } else { # No HPs
    K    = covar.fn(locations, D=D, max.dist=max.dist)
  }
  stopifnot(is(K, "Matrix"))

  expectations = EM.E(Xc, W, sigSq)

  lp   = stpca.log_posterior(Xc, K, W, rep(0,d), sigSq) # Current log posterior
  ll   = stpca.log_likelihood(Xc, W, rep(0,d), sigSq)
  dof  = d*k - 0.5*k*(k-1) + 3 + length(beta0) # Degrees of Freedom for PPCA + #HPs
  bic  = -2*ll + dof*log(n)

  stpcaObj = list(n     = nrow(Xc),
                  d     = ncol(Xc),
                  k     = k,
                  Xc    = Xc,
                  W     = W,
                  sigSq = sigSq,
                  mu    = mu,
                  V     = expectations$V,
                  Vvar  = expectations$Vvar,
                  ll    = ll,
                  lp    = lp,
                  lps   = lp,
                  bic   = bic,
                  beta  = beta0,
                  D     = D,
                  K     = K,
                  H     = list(),
                  max.dist = max.dist,
                  covar.fn = covar.fn,
                  locations = locations,
                  covar.fn.d = covar.fn.d,
                  log_evidence = NA,
                  log_evidences = c(),
                  constraints = constraints)
  class(stpcaObj) = "stpca"
  return(stpcaObj)
}
