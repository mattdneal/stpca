# All the code related to the sparsity extension lives in this file.
# The style is a little different because it is a lot newer than all
# the other R code and I've since read the Adv. R style guide.

#' Soft-thresholding operator.
#'
#' The soft-thresholding operator S_\\lambda(x) = sign(x)(abs(x) - \\lambda)_+
#'
#' @param x the vector/matrix of values to be soft thresholded
#' @param threshold the threshold lambda
soft_threshold <- function(x, threshold) {
  stopifnot(threshold>=0)
  pmax(abs(x)-threshold, 0)*sign(x)
}

w_coord_desc <- function(Xc, l, m, W, sigSq, colVmag, RtV, K, b, Kinv=NULL) {
  stopifnot(l <= ncol(W)) # l \leq k
  stopifnot(m <= nrow(W)) # m \leq d

  if (is.null(Kinv)) Kinv <- solve(K)

  # Update w_{lm}
  return(
    soft_threshold(as.numeric(
      RtV[m,l] - sigSq*Kinv[m,-m]%*%W[-m,l]
    ), sigSq/b) / (colVmag[l] + sigSq*Kinv[m,m])
  )
}

dlaplace <- function(x, b, logarithm=FALSE) {
  stopifnot(b>0)
  val <- -abs(x)/b - log(2*b)
  if (!logarithm) {val <- exp(val)}
  return(val)
}

#' The *un-normalised* prior over W, sigma^2 in SpStPCA. A product of
#' the proper, un-normalised p(W | \beta), and the improper p(sigma^2).
#'
#' @param K Prior covariance matrix
#' @param W Loadings matrix
#' @param sigSq variance of noise added to 'signal' Wv_i
#' @param b Variance on Laplace part of prior
#' @return un-normalised log prior over W
#' @export
log_sparse_prior <- function(K, W, sigSq, b) {
  return(log_sparse_prior_W(K, W, b) + log_prior_sigSq(sigSq))
}

#' The *un-normalised* prior over W in SpStPCA. This is a product of
#' the StPCA Gaussian prior and an iid Laplace prior.
#' p(W | beta, b) = \\prod^k_{i=1} N(w_i | 0, K) *
#'                  \\prod^{d}_{i=1} \\prod^{k}_{j=1} Laplace(W_ij | 0, b)
#'
#' @param K Prior covariance matrix
#' @param W Loadings matrix
#' @param b Variance on Laplace part of prior
#' @return un-normalised log prior over W
#' @export
log_sparse_prior_W <- function(K, W, b) {
  W <- Matrix(W)
  stopifnot(is.numeric(b))
  stopifnot(nrow(K) == ncol(K))
  stopifnot(nrow(W) >= ncol(W))

  stopifnot(nrow(W) >= ncol(W))
  normPart <- log_prior_W(K, W)
  laplacePart <- sum(dlaplace(W, b, logarithm=TRUE))
  return(normPart + laplacePart)
}
