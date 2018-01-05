context("EM: Expectation step")

test_that("The expectation function matches analytic solution", {
  n = 10
  d = 20
  k = 3
  Xc = scale(matrix(rnorm(n*d), nrow=n, ncol=d), scale=FALSE)
  sigSq = 1
  W     = diag(rep(1, k), nrow=d, ncol=k)

  E = EM.E(Xc, W, sigSq)

  expect_equal(E$Vmean, 0.5*Xc[,1:k])

  for (i in 1:n) {
    expect_is(E$Vvar[[i]], 'matrix')
    expect_equal(dim(E$Vvar[[i]]), c(k, k))
    expect_equal(E$Vvar[[i]], 0.5*diag(k) + 0.25*tcrossprod(Xc[i,1:k]))
  }
})

context("EM: maximisation in sigma^2")

test_that("sigma^2 maximisation stage increases complete log posterior to local max", {
  Xc     = sweep(stpca$X, 2, stpca$muHat)
  E      = EM.E(Xc, stpca$WHat, stpca$sigSqHat)
  sigSq2 = EM.M.sigSq(Xc, stpca$WHat, E$Vmean, E$Vvar)

  clp1   = complete_log_posterior(Xc, E$Vmean, E$Vvar, stpca$WHat, 0, sigSq,  stpca$K)
  clp2   = complete_log_posterior(Xc, E$Vmean, E$Vvar, stpca$WHat, 0, sigSq2, stpca$K)

  expect_gt(clp2, clp1)

  G = numDeriv::grad(function(s_) {
    complete_log_posterior(Xc, E$Vmean, E$Vvar, stpca$WHat, 0, s_, stpca$K)
  }, sigSq2)

  H = numDeriv::hessian(function(s_) {
    complete_log_posterior(Xc, E$Vmean, E$Vvar, stpca$WHat, 0, s_, stpca$K)
  }, sigSq2)

  # Zero gradient
  expect_equal(G, 0)

  # 1x1 hessian is negative definite => local maximum
  expect_lt(H, 0)
})

test_that("Iterating Expectation and sigma^2 maximisation finds likelihood local max", {
  Xc = sweep(stpca$X, 2, stpca$muHat)
  W  = stpca$WHat
  sigSq = stpca$sigSqHat

  for (i in 1:100) {
    E     = EM.E(Xc, W, sigSq)
    sigSq = EM.M.sigSq(Xc, W, E$Vmean, E$Vvar)
  }

  G = numDeriv::grad(function(sigSq_) {
    log_likelihood(Xc, W, 0, sigSq_)
  }, sigSq)[1]

  H = numDeriv::hessian(function(sigSq_) {
    log_likelihood(Xc, W, 0, sigSq_)
  }, sigSq)[1]

  # Zero gradient
  expect_equal(G, 0)

  # 1x1 hessian is negative definite => local maximum
  expect_lt(H, 0)
})

context("EM: maximisation in W")

test_that("W maximisation stage increases expected complete log posterior to zero gradient point", {
  Xc = sweep(stpca$X, 2, stpca$muHat)

  E    = EM.E(Xc, stpca$WHat, stpca$sigSqHat)

  clp1 = complete_log_posterior(Xc, E$Vmean, E$Vvar, stpca$WHat, 0,
                                stpca$sigSqHat, stpca$K)

  WHatNew = EM.M.W(Xc, stpca$sigSqHat, E$Vmean, E$Vvar, stpca$K)

  clp2 = complete_log_posterior(Xc, E$Vmean, E$Vvar, WHatNew, 0,
                                stpca$sigSqHat, stpca$K)

  expect_gte(clp2, clp1)

  #G0 = numDeriv::grad(function(W_) {
  #  complete_log_posterior(Xc, E$Vmean, E$Vvar, W_, 0, stpca$sigSqHat, stpca$K)
  #}, stpca$WHat)

  clpGrad = numDeriv::grad(function(W_) {
    complete_log_posterior(Xc, E$Vmean, E$Vvar, W_, 0, stpca$sigSqHat, stpca$K)
  }, WHatNew)

  # Zero gradient
  expect_equal(clpGrad, rep(0, length(stpca$WHat)), tol=1e-5)
})

test_that("Iterating Expectation and W maximisation finds zero gradient point", {
  require(numDeriv)
  Xc = sweep(stpca$X, 2, stpca$muHat)

  W = stpca$WHat
  for (i in 1:300) {
    E = EM.E(Xc, W, stpca$sigSqHat)
    W = EM.M.W(Xc, stpca$sigSqHat, E$Vmean, E$Vvar, stpca$K)
  }

  G = numDeriv::grad(function(W_) {
    log_likelihood(Xc, W_, 0, stpca$sigSqHat) +
    log_prior(stpca$K, W_)
  }, W)

  expect_equal(G, rep(0, length(W)), tolerance=1e-6)
})

context("EM: Full EM procedure")

test_that("fit_stpca finds the MAP theta", {
  set.seed(1)

  stpca2 = stpca$copy()
  stpca2$fit(300)

  params = list(W=stpca2$WHat, sigSq=stpca2$sigSqHat, mu=stpca2$muHat)
  veclp = function(theta) {
    L = relist(theta, params)
    lp = log_likelihood(stpca2$X, L$W, L$mu, L$sigSq) +
         log_prior(stpca2$K, L$W)
    return(lp)
  }

  thetaHat = unlist(params)
  lpMax = veclp(thetaHat)

  G = grad(veclp, thetaHat)

  expect_equal(G, rep(0, length(thetaHat)), tol=1e-7)
})
