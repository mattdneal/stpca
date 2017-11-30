context("Numerical vs analytic blocks of H")

test_that("Analytic H_{\\sigma^2} is equal to numeric H_{\\sigma^2}", {
  HsigSq.analytic = H$sigSq
  HsigSq.numeric = Matrix(numDeriv::hessian(function(sigSq_) {
    -stpca.log_posterior(X, stpcaObj$K, stpcaObj$W, stpcaObj$mu, sigSq_)
  }, x=stpcaObj$sigSq))
  expect_equal(HsigSq.analytic, HsigSq.numeric)
})

test_that("Analytic H_{\\mu} is equal to numeric H_{\\mu}", {
  Hmu.analytic = H$mu
  Hmu.numeric = Matrix(numDeriv::hessian(function(mu_) {
    -stpca.log_posterior(X, stpcaObj$K, stpcaObj$W, mu_, stpcaObj$sigSq)
  }, x=stpcaObj$mu))
  expect_equal(Hmu.analytic, Hmu.numeric, tolerance=1e-6)
})

test_that("Analytic H_{w_i} are all equal to numeric H_{w_i}", {
  for (i in 1:k) {
    Hwi.analytic = H[[paste("w", i, sep='')]]
    Hwi.numeric = Matrix(numDeriv::hessian(function(wi) {
      W_ = stpcaObj$W
      W_[,i] = wi
      -stpca.log_posterior(X, stpcaObj$K, W_, stpcaObj$mu, stpcaObj$sigSq)
    }, x=stpcaObj$W[,i]))

    expect_equal(Hwi.analytic, Hwi.numeric)
  }
})

test_that("Building H triggers an informative error if K cannot be inverted", {
  expect_error({
    with(stpcaObj, {
      K_ = cov.SE(locations, beta=log(c(1, 10000)))
      stpca.H(X, W, mu, sigSq, K_)
    })}, "Could not invert K")
})

test_that("diag of H matches numerical reconstruction on larget dataset", {
  # This is also checking non-optimised params.
  n     = 15
  k     = 4
  dim   = c(21, 21)
  d     = prod(dim)
  beta  = log(c(2, 0.2, 1e-5))
  k_se  = Curry(cov.noisy.SE, beta=beta)
  sigSq = 1.8

  synth = synthesize_data_kern(n, k, dim, kern=k_se, noisesd=sqrt(sigSq))

  Hw1diag.num = vapply(1:d, function(i) {
    as.numeric(hessian(function(wi) {
      W_ = synth$W
      W_[i,1] = wi
      -stpca.log_posterior(synth$X, synth$K, W_, 0, sigSq)
    }, synth$W[i,1]))
  }, numeric(1))

  Hw1diag = diag(stpca.H.W(synth$X, synth$W, 0, sigSq, synth$K)[[1]])

  expect_equal(Hw1diag, Hw1diag.num)
})
