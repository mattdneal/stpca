#' Calculate the log likelihood for StPCA with given parameters
#'
#' @param X Data
#' @param W
#' @param mu
#' @param sigSq
#' @return log likelihood (numeric)
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
#' l1 = log_likelihood(X, W, mu, sigSq)
#' l2 = log_likelihood(X, W%*%R, mu, sigSq)
#' stopifnot(all.equal(l1, l2))
#' @export
log_likelihood <- function(X, W, mu, sigSq) {
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

#' Compute a value proportional to the  expected complete log posterior:
#'   E[ log p(\\theta | X, V, \\beta) | V ]
#' \\propto
#'   E[ log p(X | V, \\theta) | V ] + log p(\\theta | \\beta) + E[ log p(V) | V ]
#'
#' This is "propertional to" because we do not compute the normalising
#' constant E[ log p(X, V) | V ]. This does not matter since this function is
#' used to check that each EM M-step maximizes this quantity w.r.t \\theta.
#'
#' @export
complete_log_posterior <- function(X, V, Vvar, W, mu, sigSq, K, sparse=FALSE, b=NULL) {
  require(mvtnorm)
  d = ncol(X)
  k = ncol(W)
  n = nrow(X)
  Xc = sweep(X, 2, mu)

  # E[ log p(V) | V ]
  pr1 = -0.5*(n*k*log(2*pi) +
              sum(diag(Reduce('+', Vvar))))

  # log p(\theta | \beta)
  if (sparse) {
    pr2 <- log_sparse_prior(W, K, b)
  } else {
    pr2 <- log_prior(K, W)
  }

  # E[ log p(X | V, \theta) | V ]
  pr3 = -(n*d*log(2*pi*sigSq) +
          sum(diag(crossprod(W)%*%Reduce('+', Vvar)))/sigSq -
          2*sum(vapply(1:nrow(Xc), function(n_) {
            V[n_,] %*% crossprod(W, Xc[n_,])
          }, numeric(1)))/sigSq +
          (norm(Xc, 'F')^2)/sigSq)*0.5

  return(pr1 + pr2 + pr3)
}
