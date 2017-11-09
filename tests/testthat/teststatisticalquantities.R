context("statistical quantities")

test_that("The posterior is unaffected by adding a vector to colMeans(X) and mu", {
  lp1 = with(stpcaObj, stpca.log_posterior(X, K, W, mu, sigSq))
  lp2 = with(stpcaObj, stpca.log_posterior(Xc, K, W, rep(0,ncol(Xc)), sigSq))
  expect_equal(lp1, lp2)
})

test_that("The posterior is unaffected by rotating W", {
  set.seed(1)
  lp1 = with(stpcaObj, stpca.log_posterior(X, K, W, mu, sigSq))

  for (i in 1:10) {
    R = svd(matrix(rnorm(k*k), ncol=k))$u # random orthonornal matrix
    lp2 = with(stpcaObj, stpca.log_posterior(Xc, K, W%*%R, rep(0,ncol(Xc)), sigSq))
    expect_equal(lp1, lp2)
  }
})

test_that("The maximum of the likelihood is given by PPCA", {
  ## PPCA MAP solution
  mu    = colMeans(X)
  Xc    = scale(X, scale=FALSE, center=TRUE)
  X.eig = svd(Xc/sqrt(n), nu=0, nv=k)
  sigSq = sum(X.eig$d[-(1:k)]^2)/(d-k)
  W     = X.eig$v %*% diag(sqrt(X.eig$d[1:k]^2 - sigSq), nrow=k, ncol=k)

  params = list(W=W, sigSq=sigSq, mu=mu)
  vecll = function(theta) {
    L = relist(theta, params)
    return(stpca.log_likelihood(X, L$W, L$mu, L$sigSq))
  }

  thetaHat = unlist(params)
  llMax = vecll(thetaHat)

  # All perturbations to the maximum likelihood parameters cause a lower log
  # likelihood
  for (i in 1:50) {
    llPerturb = vecll(unlist(params) + rnorm(length(thetaHat))*1e-4)
    expect_true(llMax > llPerturb)
  }
})

test_that("Prior matches simple analytical solution with W=0", {
  # W=0 means log prior is just 0.5*(log det 2*pi*K )
  W = matrix(0, nrow=d, ncol=k)
  lpr = stpca.log_prior(stpcaObj$K, W)
  logDetK = as.numeric(determinant(stpcaObj$K, log=TRUE)$modulus)
  lpr.analytic = -0.5*k*(d*log(2*pi) + logDetK)
  expect_equal(lpr, lpr.analytic)
})

test_that("Prior matches simple analytical solution with K=I (dense)", {
  # K = I means prior is just an independand normal on each element of W.
  K = diag(d)
  for (i in 1:10) {
    W = matrix(rnorm(d*k), nrow=d)
    lpr = stpca.log_prior(K, W)
    lpr.analytic = -(d*k*log(2*pi) + sum(W^2))/2
    expect_equal(lpr, lpr.analytic)
  }
})

test_that("Prior matches simple analytical solution with K=I (sparse)", {
  K = Matrix(diag(d)) # Test sparsity
  for (i in 1:10) {
    W = matrix(rnorm(d*k), nrow=d)
    lpr = stpca.log_prior(K, W)
    lpr.analytic = -(d*k*log(2*pi) + sum(W^2))/2
    expect_equal(lpr, lpr.analytic)
  }
})

test_that("Likelihood matches simple analytical solution with W=0", {
  # W=0 means log likelihood
  W = matrix(0, nrow=d, ncol=k)
  mu = rep(0, d)
  sigSq = 1
  ll = stpca.log_likelihood(X, W, mu, sigSq)

  ll.analytic = sum(dnorm(c(X), mean=mu, sd=sigSq, log=TRUE))

  expect_equal(ll, ll.analytic)
})
