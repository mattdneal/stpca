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
log_evidence <- function(X, K, W, mu, sigSq) {
  if (!any(is.finite(diag(K)))) {
    return(-Inf)
  }

  n = nrow(X)
  d = ncol(X)
  k = ncol(W)

  # If the inversion cannot be done, logZ defaults to -Inf
  logZ = -Inf
  try({
    H = compute_H(X, W, mu, sigSq, K)
    logDetH = sum(vapply(H, function(Hblock) {
       as.numeric(determinant(Hblock, logarithm=TRUE)$modulus)
    }, numeric(1)))

    # Laplace-approximated log evidence
    logZ = (log_prior(K, W) +
            log_likelihood(X, W, mu, sigSq) +
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
log_evidence_d <- function(X, K, W, mu, sigSq, beta, dK) {
  success=FALSE
  try({HW = H.W(X, W, mu, sigSq, K); success=TRUE})

  if(!success) {return(lapply(seq_along(beta), function(b) 0))}

  logPrior.d = log_prior_d(W, beta, K, dK)
  logDetH.d  = log_det_H_d(K, dK, HW)

  logEvidence.d = vapply(seq_along(beta), function(i) {
    logPrior.d[i] - 0.5*logDetH.d[i]
  }, numeric(1))
  names(logEvidence.d) = names(beta)
  return(logEvidence.d)
}
