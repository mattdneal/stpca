context("Numerical vs analytic blocks of H")

test_that("Analytic H_{\\sigma^2} is equal to numeric H_{\\sigma^2}", {
  HsigSq.analytic = H$sigSq
  HsigSq.numeric = Matrix(numDeriv::hessian(function(sigSq_) {
    -log_likelihood(X, stpca$WHat, stpca$muHat, sigSq_)
  }, x=stpca$sigSqHat))
  expect_equal(HsigSq.analytic, HsigSq.numeric)
})

test_that("Analytic H_{\\mu} is equal to numeric H_{\\mu}", {
  Hmu.analytic = H$mu
  Hmu.numeric = Matrix(numDeriv::hessian(function(mu_) {
    -(log_likelihood(stpca$X, stpca$WHat, mu_, stpca$sigSqHat) +
      log_prior(stpca$K, stpca$WHat))
  }, x=stpca$muHat))
  expect_equal(Hmu.analytic, Hmu.numeric, tolerance=1e-6)
})

test_that("Analytic H_{w_i} are all equal to numeric H_{w_i}", {
  for (i in 1:k) {
    Hwi.analytic = H[[paste("w", i, sep='')]]
    Hwi.numeric = numDeriv::hessian(function(wi) {
      W_ = stpca$WHat
      W_[,i] = wi
      -(log_likelihood(stpca$X, W_, stpca$muHat, stpca$sigSqHat) +
        log_prior(stpca$K, W_))
    }, x=stpca$WHat[,i])

    expect_equal(as.matrix(Hwi.analytic), Hwi.numeric, tolerance=1e-6,
                 check.attributes=FALSE)
  }
})

test_that("Building H triggers an informative error if K cannot be inverted", {
  expect_error({
    K_ = cov.SE(locs, beta=log(c(1, 10000)))
    compute_H(X, stpca$WHat, stpca$muHat, stpca$sigSqHat, K_)
  }, "Could not invert K")
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
      -(log_likelihood(synth$X, W_, 0, sigSq) +
        log_prior(synth$K, W_))
    }, synth$W[i,1]))
  }, numeric(1))

  Hw1diag = diag(compute_H_W(synth$X, synth$W, 0, sigSq, synth$K)[[1]])

  expect_equal(Hw1diag, Hw1diag.num)
})
