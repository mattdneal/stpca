#' De-noise a sample using a trained \code{stpca} object.
#'
#' @param object An \code{stpca} object returned from a call to \code{stpca}
#' @param samples samples to de-noise
#' @return De-noised samples of the same dimensionality as the parameter \code{samples}
#' @export
predict.stpca <- function(object, ...) {
  samples = list(...)[[1]]

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
