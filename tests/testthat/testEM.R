context("expectation")

test_that("The expectation function matches analytic solution", {
  n = 10
  d = 20
  k = 3
  Xc = scale(matrix(rnorm(n*d), nrow=n, ncol=d), scale=FALSE)
  sigSq = 1
  W     = diag(rep(1, k), nrow=d, ncol=k)

  expectations = EM.E(Xc, W, sigSq)

  expect_equal(expectations$V, 0.5*Xc[,1:k])

  for (i in 1:n) {
    expect_is(expectations$Vvar[[i]], 'matrix')
    expect_equal(dim(expectations$Vvar[[i]]), c(k, k))
    expect_equal(expectations$Vvar[[i]], 0.5*diag(k) + 0.25*tcrossprod(Xc[i,1:k]))
  }
})

## Generate the data
library(functional)

n     = 15
k     = 4
dim   = c(11, 11)
d     = prod(dim)
beta  = log(c(2, 0.2))
k_se  = Curry(cov.SE, beta=beta)
sigSq = 1.8

dataset = synthesize_data_kern(n, k, dim, kern=k_se, noisesd=sqrt(sigSq))
locations = dataset$grid
X = dataset$X

beta0 = cov.SE.beta0(X, locations, k)

stpcaObj = stpca(X, k, locations, cov.SE, cov.SE.d, beta0, trace=0, maxit.inner=0, maxit.outer=0)
#stpcaObj = stpca.iterate.beta(stpcaObj) # Update beta, now the thetas need updating

context("maximisation in sigma^2")

test_that("sigma^2 maximisation stage increases complete log posterior", {
  set.seed(1)
  stpcaObj$W = matrix(rnorm(d*k), nrow=d, ncol=k)
  clp1     = with(stpcaObj, stpca.complete_log_posterior(X, V, W, mu, sigSq, K))
  E        = with(stpcaObj, EM.E(Xc, W, sigSq))
  sigSqNew = with(stpcaObj, EM.M.sigSq(Xc, W, E$V, E$Vvar))
  clp2     = with(stpcaObj, stpca.complete_log_posterior(X, V, W, mu, sigSqNew, K))
  expect_gt(clp2, clp1)
})

test_that("Iterating Expectation and sigma^2 maximisation finds local max", {
  Xc = stpcaObj$Xc
  K = stpcaObj$K
  W = stpcaObj$W
  mu = rep(0, d)
  sigSq = stpcaObj$sigSq
  for (i in 1:200) {
    E     = EM.E(Xc, stpcaObj$W, sigSq)
    sigSq = EM.M.sigSq(Xc, stpcaObj$W, E$V, E$Vvar)
  }

  lp = stpca.log_posterior(Xc, K, W, mu, sigSq)

  delta = 1e-8
  lpUp = stpca.log_posterior(Xc, K, W, mu, sigSq+delta)
  lpDn = stpca.log_posterior(Xc, K, W, mu, sigSq-delta)

  #gradient = (lpUp-lpDn)/(2*delta)
  #expect_lt(gradient, 1e-6)

  expect_gte(lp, lpUp)
  expect_gte(lp, lpDn)
})

context("maximisation in W")

test_that("W maximisation stage increases complete log posterior", {
  clp1 = with(stpcaObj, stpca.complete_log_posterior(Xc, V, W, rep(0,d), sigSq, K))
  Wnew  = with(stpcaObj, EM.M.W(Xc, sigSq, V, Vvar, K))
  clp2 = with(stpcaObj, stpca.complete_log_posterior(Xc, V, Wnew, rep(0,d), sigSq, K))
  expect_gte(clp2, clp1)
})

test_that("Iterating Expectation and W maximisation finds local max", {
  require(numDeriv)
  Xc = stpcaObj$Xc
  K = stpcaObj$K
  W = stpcaObj$W
  mu = rep(0, d)
  sigSq = stpcaObj$sigSq
  for (i in 1:200) {
    E = EM.E(Xc, W, sigSq)
    W = EM.M.W(Xc, sigSq, E$V, E$Vvar, K)
  }

  lp = stpca.log_posterior(Xc, K, W, mu, sigSq)

  for (i in 1:300) {
    Wpert = W + 1e-5*matrix(rnorm(nrow(W)*ncol(W)), nrow=nrow(W))
    lpPert = stpca.log_posterior(Xc, K, W, mu, sigSq)
    expect_gte(lp, lpPert)
  }
})
