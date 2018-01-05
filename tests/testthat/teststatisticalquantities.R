context("statistical quantities")

test_that("The likelihood is unaffected by adding a vector to colMeans(X) and mu", {
  ll1 = log_likelihood(stpca$X, stpca$WHat, stpca$muHat, stpca$sigSqHat)
  ll2 = log_likelihood(scale(stpca$X, scale=FALSE),
                       stpca$WHat, 0, stpca$sigSqHat)
  expect_equal(ll1, ll2)
})

test_that("The log likelihood is unaffected by rotating W", {
  set.seed(1)

  ll1 = log_likelihood(stpca$X, stpca$WHat, stpca$muHat, stpca$sigSqHat)

  for (i in 1:10) {
    R = svd(matrix(rnorm(k*k), ncol=k))$u # random orthonornal matrix
    ll2 = log_likelihood(stpca$X, stpca$WHat%*%R, stpca$muHat, stpca$sigSqHat)
    expect_equal(ll1, ll2)
  }
})

test_that("The maximum of the likelihood is given by PPCA", {
  ## PPCA MAP solution
  mu    = colMeans(stpca$X)
  Xc    = scale(stpca$X, scale=FALSE, center=TRUE)
  X.eig = svd(Xc/sqrt(stpca$n), nu=0, nv=k)
  sigSq = sum(X.eig$d[-(1:stpca$k)]^2)/(stpca$d-stpca$k)
  W     = X.eig$v %*% diag(sqrt(X.eig$d[seq_len(stpca$k)]^2 - sigSq))

  params = list(W=W, sigSq=sigSq, mu=mu)
  vecll = function(theta) {
    L = relist(theta, params)
    return(log_likelihood(stpca$X, L$W, L$mu, L$sigSq))
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
  lpr = log_prior(stpca$K, W)
  logDetK = as.numeric(determinant(stpca$K, log=TRUE)$modulus)
  lpr.analytic = -0.5*k*(d*log(2*pi) + logDetK)
  expect_equal(lpr, lpr.analytic)
})

test_that("Prior matches simple analytical solution with K=I (dense)", {
  # K = I means prior is just an independand normal on each element of W.
  K = diag(d)
  for (i in 1:10) {
    W = matrix(rnorm(d*k), nrow=d)
    lpr = log_prior(K, W)
    lpr.analytic = -(d*k*log(2*pi) + sum(W^2))/2
    expect_equal(lpr, lpr.analytic)
  }
})

test_that("Prior matches simple analytical solution with K=I (sparse)", {
  K = Matrix(diag(d)) # Test sparsity
  for (i in 1:10) {
    W = matrix(rnorm(d*k), nrow=d)
    lpr = log_prior(K, W)
    lpr.analytic = -(d*k*log(2*pi) + sum(W^2))/2
    expect_equal(lpr, lpr.analytic)
  }
})

test_that("Log prior matches a multivariate normal", {
  require(mvtnorm)
  lp1 = log_prior(stpca$K, stpca$WHat)
  lp2 = sum(dmvnorm(t(stpca$WHat), sigma=as.matrix(stpca$K), log=TRUE))
  expect_equal(lp1, lp2)
})

test_that("Likelihood matches simple analytical solution with W=0", {
  # W=0 means log likelihood
  W = matrix(0, nrow=d, ncol=k)
  mu = rep(0, d)
  sigSq = 1
  ll = log_likelihood(stpca$X, W, mu, sigSq)

  ll.analytic = sum(dnorm(c(X), mean=mu, sd=sigSq, log=TRUE))

  expect_equal(ll, ll.analytic)
})
