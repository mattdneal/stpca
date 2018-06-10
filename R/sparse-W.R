# All the code related to the sparsity extension lives in this file.
# The style is a little different because it is a lot newer than all
# the other R code and I've since read the Adv. R style guide.

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

log_sparse_prior <- function(W, K, b) {
  W <- Matrix(W)
  stopifnot(nrow(W) >= ncol(W))
  normPart <- log_prior(K, W)
  laplacePart <- sum(dlaplace(W, b, logarithm=TRUE))
  return(normPart + laplacePart)
}
