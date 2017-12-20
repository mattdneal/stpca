#' Synthesize fake data from StPCA model
#' @export
#' @param n Number of samples to create
#' @param k Latent dimensionality
#' @param K prior covariance matrix
#' @param noisesd standard deviation of noise added to data
#' @import MASS
synthesize_data <- function(n, k, K, noisesd=0) {
  d   = ncol(K)
  stopifnot(k<=d)
  R_K = chol(K)
  W   = t(R_K) %*% matrix(rnorm(d*k), nrow=d, ncol=k)
  V   = mvrnorm(n=n, mu=rep(0, k), Sigma=diag(k))
  X   = V %*% t(W) + matrix(rnorm(n*d, sd=noisesd), nrow=n, ncol=d)

  W.svd = svd(W)
  W     = W %*% W.svd$v
  V     = V %*% W.svd$v

  return(list(X   = X,
              W   = W,
              V   = V))
}

#' Synthesize fake data from StPCA model
#' @export
#' @param n Number of samples to create
#' @param k Latent dimensionality
#' @param dim dimensions of grid to construct
#' @param kern cvarinace function used to build the prior covariance matrix K
#' @param noisesd standard deviation of noise added to data
#' @import fields
synthesize_data_kern <- function(n, k, dim, kern, noisesd=0) {
  grid = as.matrix(make.surface.grid(lapply(dim, function(d_) seq(0, 1, length=d_))))
  K    = kern(grid)
  data.synth = synthesize_data(n, k, K, noisesd)
  return(list(X    = data.synth$X,
              W    = data.synth$W,
              V    = data.synth$V,
              grid = grid,
              K    = K))
}
