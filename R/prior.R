#' Calculate the *un-normalised* log prior (only in sigSq) for StPCA
#' with the given loadings matrix
#'
#' @param K Prior covariance matrix
#' @param W Loadings matrix
#' @return un-normalised log prior (numeric)
#' @export
log_prior <- function(K, W, sigSq) {
  return(log_prior_W(K, W) + log_prior_sigSq(sigSq))
}

#' The improper prior over sigSq. Proportional to sigma^-2
#'
#' @param sigSq
#' @return un-normalised log prior (numeric)
#' @export
log_prior_sigSq <- function(sigSq) {
  stopifnot(is.finite(sigSq))
  stopifnot(sigSq>0)
  return(-log(sigSq))
}

#' The proper prior over W.
#' p(W) = \\prod^k_{i=1} N(w_i | 0, K)
#'
#' @param K Prior covariance matrix
#' @param W Loadings matrix
#' @return normalised log prior (numeric)
#' @export
log_prior_W <- function(K, W) {
  W <- Matrix(W)
  stopifnot(nrow(K) == ncol(K))
  stopifnot(nrow(W) >= ncol(W))

  d = nrow(W)
  k = ncol(W)

  # This special case for sparse matrices is more numerically stable
  if (is(K, "sparseMatrix")) {
    if (any(!is.finite(diag(K)))) { return(-Inf) }

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
    if (any(!is.finite(diag(K)))) { return(-Inf) }
    trWtKinvW = sum(diag(Matrix::crossprod(W, solve(K, W)))) # Fnorm(R\W)^2?
    # TODO: Calculate determinant from decomposition
    logDetK = as.numeric(determinant(K, logarithm=TRUE)$modulus)
  }

  return(-0.5*( k*d*log(2*pi) +
                k*logDetK +
                trWtKinvW))
}

#' Compute the partial derivatives of the log prior with respect to the
#' hyperparameters with the value of beta provided.
#'
#' @param W Loadings matrix
#' @param beta hyperparameter values
#' @param K Prior covariance matrix
#' @param KD Prior covariance matrix derivatives
#' @return Partial derivatives of log prior
#' @export
log_prior_d <- function(W, beta, K, KD) {
  deriv = numeric(length(beta))
  KinvW = solve(K, W)

  for (i in seq_along(beta)) {
    term1 = ncol(W) * sum(diag( solve(K, KD[[i]])  ))

    term2 = -sum(diag( tcrossprod(KinvW)%*%KD[[i]] ))
    #term2b = -sum(vapply(1:ncol(W), function(k_) {
    #  as.numeric(KinvW[,k_] %*% KD[[i]] %*% KinvW[,k_])
    #}, numeric(1)))

    deriv[i] = -0.5*(term1 + term2)
  }
  names(deriv) = names(beta)
  return(deriv)
}
