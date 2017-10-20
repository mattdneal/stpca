#' Calculate the log likelihood for StPCA with given parameters
#'
#' @param X Data
#' @param W
#' @param mu
#' @param sigSq
#' @return log likelihood (numeric)
#' @export
#' @examples
#' d = 50
#' k = 5
#' n = 15
#'
#' set.seed(1)
#' X  = matrix(rnorm(n*d), nrow=n, ncol=d)
#' mu = colMeans(X)
#' Xc = sweep(X, 2, mu, '-')
#' covar.svd = svd(Xc/sqrt(n), nu=0, nv=k)
#' covar.eigval = covar.svd$d^2
#' sigSq = sum(covar.eigval[-(1:k)])/(d-k)
#' W     = covar.svd$v %*% diag(sqrt(covar.eigval[1:k] - sigSq), ncol=k, nrow=k)
#'
#' R     = svd(matrix(rnorm(k*k), ncol=k, nrow=k))$u # Random orthonormal matrix
#'
#' # The likelihood is invariant to multiplying  by an orthonormal matrix.
#' l1 = stpca.log_likelihood(X, W, mu, sigSq)
#' l2 = stpca.log_likelihood(X, W%*%R, mu, sigSq)
#' stopifnot(all.equal(l1, l2))
stpca.log_likelihood <- function(X, W, mu, sigSq) {
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

#' Calculate the *un-normalised* log prior for StPCA with given loadings matrix
#'
#' @param K Prior covariance matrix
#' @param W Loadings matrix
#' @return un-normalised log prior (numeric)
#' @export
stpca.log_prior <- function(K, W) {
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

#' Calculate the *un-normalised* log posterior for StPCA with given parameters
#'
#' @param X Data
#' @param K Prior covariance matrix
#' @param W Loadings matrix
#' @param mu
#' @param sigSq
#' @return un-normalised log posterior (numeric)
#' @export
stpca.log_posterior <- function(X, K, W, mu, sigSq) {
  return(stpca.log_likelihood(X, W, mu, sigSq) + stpca.log_prior(K, W))
}
