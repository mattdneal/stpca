#' Initialize mu, sigSq and W from PPCA.
#'
#' @param X Data
#' @param k latent dimensionality
#'
#' @return MAP PPCA parameters
initialize_from_ppca <- function(X, k) {
  n <- nrow(X)
  d <- ncol(X)

  Xc <- scale(X, scale=FALSE) # Centered data: nxd
  covar.svd <- svd(Xc/sqrt(n), nu=0, nv=k)
  covar.eigval <- covar.svd$d^2


  mu    = attr(Xc, "scaled:center")
  sigSq = sum(covar.eigval[-(1:k)])/(d-k)
  W     = (covar.svd$v %*%
           diag(sqrt(covar.eigval[seq_len(k)] - sigSq), nrow=k, ncol=k))

  return(list(mu=mu, sigSq=sigSq, W=W))
}
