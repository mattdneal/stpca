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

context("maximisation in sigma^2")

test_that("sigma^2 maximisation stage increases complete log posterior to local max", {
  E      = with(stpcaObj, EM.E(Xc, W, sigSq))
  sigSq2 = with(stpcaObj, EM.M.sigSq(Xc, W, E$V, E$Vvar))
  clp1   = with(stpcaObj, stpca.complete_log_posterior(Xc, E$V, E$Vvar, W, 0, sigSq,  K))
  clp2   = with(stpcaObj, stpca.complete_log_posterior(Xc, E$V, E$Vvar, W, 0, sigSq2, K))

  expect_gt(clp2, clp1)

  G = numDeriv::grad(function(s_) {
    with(stpcaObj, stpca.complete_log_posterior(Xc, E$V, E$Vvar, W, 0, s_, K))
  }, sigSq2)

  H = numDeriv::hessian(function(s_) {
    with(stpcaObj, stpca.complete_log_posterior(Xc, E$V, E$Vvar, W, 0, s_, K))
  }, sigSq2)

  # Zero gradient
  expect_equal(G, 0)

  # 1x1 hessian is negative definite => local maximum
  expect_lt(H, 0)
})

test_that("Iterating Expectation and sigma^2 maximisation finds posterior local max", {
  Xc = stpcaObj$Xc
  K = stpcaObj$K
  W = stpcaObj$W
  sigSq = stpcaObj$sigSq
  for (i in 1:200) {
    E     = EM.E(Xc, W, sigSq)
    sigSq = EM.M.sigSq(Xc, stpcaObj$W, E$V, E$Vvar)
  }

  G = numDeriv::grad(function(sigSq_) {
    stpca.log_posterior(Xc, K, W, 0, sigSq_)
  }, sigSq)[1]

  H = numDeriv::hessian(function(sigSq_) {
    stpca.log_posterior(Xc, K, W, 0, sigSq_)
  }, sigSq)[1]

  # Zero gradient
  expect_equal(G, 0)

  # 1x1 hessian is negative definite => local maximum
  expect_lt(H, 0)
})

context("maximisation in W")

test_that("W maximisation stage increases complete log posterior to zero gradient point", {
  E    = with(stpcaObj, EM.E(Xc, W, sigSq))
  clp1 = with(stpcaObj, stpca.complete_log_posterior(Xc, E$V, E$Vvar, W, rep(0,d), sigSq, K))
  Wnew = with(stpcaObj, EM.M.W(Xc, sigSq, E$V, E$Vvar, K))
  clp2 = with(stpcaObj, stpca.complete_log_posterior(Xc, E$V, E$Vvar, Wnew, rep(0,d), sigSq, K))
  expect_gte(clp2, clp1)

  G = numDeriv::grad(function(W_) {
    with(stpcaObj, stpca.complete_log_posterior(Xc, E$V, E$Vvar, W_, 0, sigSq, K))
  }, Wnew)

  # Zero gradient
  expect_equal(G, rep(0, length(stpcaObj$W)))
})

test_that("Iterating Expectation and W maximisation finds zero gradient point", {
  require(numDeriv)
  Xc = stpcaObj$Xc
  K = stpcaObj$K
  W = stpcaObj$W
  mu = rep(0, d)
  sigSq = stpcaObj$sigSq
  for (i in 1:300) {
    E = EM.E(Xc, W, sigSq)
    W = EM.M.W(Xc, sigSq, E$V, E$Vvar, K)
  }

  G = numDeriv::grad(function(W_) {
    stpca.log_posterior(Xc, K, W_, 0, sigSq)
  }, W)

  expect_equal(G, rep(0, length(W)), tolerance=1e-6)

  lp = stpca.log_posterior(Xc, K, W, mu, sigSq)
  for (i in 1:300) {
    Wpert = W + 1e-5*matrix(rnorm(nrow(W)*ncol(W)), nrow=nrow(W))
    lpPert = stpca.log_posterior(Xc, K, Wpert, mu, sigSq)
    expect_gte(lp, lpPert)
  }
})

context("Full EM procedure")

test_that("stpca.iterate.theta maximises the log posterior", {
  set.seed(1)
  stpcaObj = stpca.iterate.theta(stpcaObj, maxit.inner=50)
  lp = with(stpcaObj, stpca.log_posterior(Xc, K, W, rep(0,ncol(Xc)), sigSq))

  for (i in 1:800) {
    Wpert  = stpcaObj$W + 1e-5*matrix(rnorm(d*k), nrow=d, ncol=k)
    sspert = stpcaObj$sigSq + 1e-5*rnorm(1)
    lpPert = with(stpcaObj, stpca.log_posterior(Xc, K, Wpert, mu, sspert))
    expect_gt(lp, lpPert)
  }
})
