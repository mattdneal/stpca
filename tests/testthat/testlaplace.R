context("Laplace approximation")

set.seed(1)

library(functional)

n     = 15
k     = 2
dim   = c(5, 5)
d     = prod(dim)
beta.real = log(c(2, 0.2))
k_se  = Curry(cov.SE, beta=beta.real)
sigSq = 1.8

dataset = synthesize_data_kern(n, k, dim, kern=k_se, noisesd=sqrt(sigSq))
locations = dataset$grid
X = dataset$X

beta0 = cov.SE.beta0(X, locations, k)

## Fit the data
stpcaObj = stpca.init(X, k, locations, cov.SE, cov.SE.d, beta0, trace=0)

H = stpca.H(X, stpcaObj$W, stpcaObj$mu, stpcaObj$sigSq, stpcaObj$K)

test_that("Analytic H_{\\sigma^2} is equal to numeric H_{\\sigma^2}", {
  HsigSq.analytic = H$sigSq
  HsigSq.numeric = Matrix(hessian(function(sigSq_) {
    -stpca.log_posterior(X, stpcaObj$K, stpcaObj$W, stpcaObj$mu, sigSq_)
  }, x=stpcaObj$sigSq))
  expect_equal(HsigSq.analytic, HsigSq.numeric)
})

test_that("Analytic H_{\\mu} is equal to numeric H_{\\mu}", {
  Hmu.analytic = H$mu
  Hmu.numeric = Matrix(hessian(function(mu_) {
    -stpca.log_posterior(X, stpcaObj$K, stpcaObj$W, mu_, stpcaObj$sigSq)
  }, x=stpcaObj$mu))
  expect_equal(Hmu.analytic, Hmu.numeric, tolerance=1e-6)
})

test_that("Analytic H_{w_i} are all equal to numeric H_{w_i}", {
  for (i in 1:k) {
    Hwi.analytic = H[[paste("w", i, sep='')]]
    Hwi.numeric = Matrix(hessian(function(wi) {
      W = stpcaObj$W
      W[,i] = wi
      -stpca.log_posterior(X, stpcaObj$K, W, stpcaObj$mu, stpcaObj$sigSq)
    }, x=stpcaObj$W[,i]))

    expect_equal(Hwi.analytic, Hwi.numeric, tolerance=1e-6)
  }
})
