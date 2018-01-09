#' Compute the laplace approximation to the log evidence given the MAP
#' parameters K, mu, sigSq as well as the prior covariance matrix K.
#' Note that this is multiplied by an UN-KNOWN CONSTANT due to the flat
#' priors over mu and sigSq. However, this unknown constant is always
#' the same regardless of k and K, so this may be used to compute
#' meaningful bayes factors between StPCA models.
#'
#' @param X Data
#' @param K Prior covariance matrix
#' @param WHat Loadings matrix
#' @param muHat
#' @param sigSqHat
#' @return Approximate log evidence
log_evidence <- function(X, K, WHat, muHat, sigSqHat, H) {
  if (!any(is.finite(diag(K)))) {
    return(-Inf)
  }

  n = nrow(X)
  d = ncol(X)
  k = ncol(WHat)

  logPrior = log_prior(K, WHat)

  logLik = log_likelihood(X, WHat, muHat, sigSqHat)

  logDetH = sum(vapply(H, function(Hblock) {
     as.numeric(determinant(Hblock, logarithm=TRUE)$modulus)
  }, numeric(1)))

  # Laplace-approximated log evidence
  logZ = (logPrior +
          logLik +
          (0.5*(d*k+d+1))*log(2*pi) -
          0.5*logDetH)
  return(logZ)
}

#' Compute the derivative of the approximate log evidence with respect to the
#' value of beta provided.
#'
#' @param X Data
#' @param K Prior covariance matrix
#' @param WHat Loadings matrix
#' @param muHat
#' @param sigSqHat
#' @param beta
#' @param KD
#' @return Partial derivatives of approximate log evidence
log_evidence_d <- function(X, K, WHat, muHat, sigSqHat, beta, KD, H) {
  #HW = H[grep("^w", names(H))]

  success=FALSE
  try({
    HW = compute_H_W(X, WHat, muHat, sigSqHat, K)
    success=TRUE
  })

  logPriorD = log_prior_d(WHat, beta, K, KD)
  logDetHD  = log_det_H_d(K, KD, HW)

  logEvidenceD = vapply(seq_along(beta), function(i) {
    logPriorD[i] - 0.5*logDetHD[i]
  }, numeric(1))
  names(logEvidenceD) = names(beta)
  return(logEvidenceD)
}
