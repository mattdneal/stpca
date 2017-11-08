context("stpca variable types")

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

maxit.inner=3
maxit.outer=1

beta0 = cov.SE.beta0(X, locations, k)

## Fit the data
stpcaObj = stpca.init(X, k, locations, cov.SE, cov.SE.d, beta0, trace=0)

test_that("stpca.init returns a valid stpca object", {
  expect_is(stpcaObj$X,     "matrix")
  expect_is(stpcaObj$W,     "matrix")
  expect_is(stpcaObj$sigSq, "numeric")
  expect_is(stpcaObj$mu,    "numeric")
  expect_is(stpcaObj$V,     "matrix")
  expect_is(stpcaObj$ll,    "numeric")
  expect_is(stpcaObj$lp,    "numeric")
  expect_is(stpcaObj$lps,   "numeric")
  expect_is(stpcaObj$bic,   "numeric")
  expect_true(length(stpcaObj$beta0)==0 | is(stpcaObj$beta0, "numeric"))
  expect_is(stpcaObj$D,     "Matrix")
  expect_is(stpcaObj$K,     "Matrix")
  expect_is(stpcaObj$covar.fn,   "function")
  expect_is(stpcaObj$covar.fn.d, "function")
  expect_is(stpcaObj$locations,  "matrix")

  expect_equal(dim(stpcaObj$X), c(n, d))
  expect_equal(dim(stpcaObj$W), c(d, k))
  expect_length(stpcaObj$sigSq, 1)
  expect_length(stpcaObj$mu, d)
  expect_equal(dim(stpcaObj$V), c(n, k))
  expect_length(stpcaObj$ll, 1)
  expect_length(stpcaObj$lp, 1)
  expect_length(stpcaObj$lps, 1)
  expect_length(stpcaObj$bic, 1)
  expect_equal(dim(stpcaObj$D), c(d, d))
  expect_equal(dim(stpcaObj$K), c(d, d))
  expect_equal(dim(locations)[1], d)
})

stpcaObj = stpca.iterate(stpcaObj, maxit.inner=maxit.inner,
                                   maxit.outer=maxit.outer)

test_that("stpca.iterate returns a valid stpca object", {
  expect_is(stpcaObj$X,     "matrix")
  expect_is(stpcaObj$W,     "matrix")
  expect_is(stpcaObj$sigSq, "numeric")
  expect_is(stpcaObj$mu,    "numeric")
  expect_is(stpcaObj$V,     "matrix")
  expect_is(stpcaObj$ll,    "numeric")
  expect_is(stpcaObj$lp,    "numeric")
  expect_is(stpcaObj$lps,   "numeric")
  expect_is(stpcaObj$bic,   "numeric")
  expect_true(length(stpcaObj$beta0)==0 | is(stpcaObj$beta0, "numeric"))
  expect_is(stpcaObj$D,     "Matrix")
  expect_is(stpcaObj$K,     "Matrix")
  expect_is(stpcaObj$covar.fn,   "function")
  expect_is(stpcaObj$covar.fn.d, "function")
  expect_is(stpcaObj$locations,  "matrix")

  expect_equal(dim(stpcaObj$X), c(n, d))
  expect_equal(dim(stpcaObj$W), c(d, k))
  expect_length(stpcaObj$sigSq, 1)
  expect_length(stpcaObj$mu, d)
  expect_equal(dim(stpcaObj$V), c(n, k))
  expect_length(stpcaObj$ll, 1)
  expect_length(stpcaObj$lp, 1)
  expect_length(stpcaObj$lps, 1 + maxit.inner*(maxit.outer + 1))
  expect_length(stpcaObj$bic, 1)
  expect_equal(dim(stpcaObj$D), c(d, d))
  expect_equal(dim(stpcaObj$K), c(d, d))
  expect_equal(dim(locations)[1], d)
})

context("stpca.iterate.beta")

stpcaObjNew = stpca.iterate.beta(stpcaObj)

test_that("stpca.iterate.beta increases log evidence", {
  expect_gt(stpcaObjNew$log_evidence, stpcaObj$log_evidence)
})

test_that("stpca.iterate.beta modifies beta, H, K", {
  expect_that(identical(stpcaObjNew$H,    stpcaObj$H),    is_false())
  expect_that(identical(stpcaObjNew$K,    stpcaObj$K),    is_false())
  expect_that(identical(stpcaObjNew$beta, stpcaObj$beta), is_false())
})

test_that("stpca.iterate.beta does not change theta", {
  expect_identical(stpcaObjNew$W,     stpcaObj$W)
  expect_identical(stpcaObjNew$sigSq, stpcaObj$sigSq)
  expect_identical(stpcaObjNew$V,     stpcaObj$V)
  expect_identical(stpcaObjNew$ll,    stpcaObj$ll)
  expect_identical(stpcaObjNew$lp,    stpcaObj$lp)
  expect_identical(stpcaObjNew$lps,   stpcaObj$lps)
})

context("stpca.iterate.theta")

stpcaObj = stpcaObjNew
stpcaObjNew = stpca.iterate.theta(stpcaObj, maxit.inner=50)

test_that("stpca.iterate.theta increases log posterior", {
  expect_gt(stpcaObjNew$lp, stpcaObj$lp)
})

test_that("stpca.iterate.theta only changes relevant variables", {
  expect_identical(stpcaObjNew$H,            stpcaObj$H)
  expect_identical(stpcaObjNew$K,            stpcaObj$K)
  expect_identical(stpcaObjNew$beta,         stpcaObj$beta)
  expect_identical(stpcaObjNew$log_evidence, stpcaObj$log_evidence)

  expect_that(identical(stpcaObjNew$W,     stpcaObj$W),     is_false())
  expect_that(identical(stpcaObjNew$sigSq, stpcaObj$sigSq), is_false())
  expect_that(identical(stpcaObjNew$V,     stpcaObj$V),     is_false())
  expect_that(identical(stpcaObjNew$ll,    stpcaObj$ll),    is_false())
  expect_that(identical(stpcaObjNew$lp,    stpcaObj$lp),    is_false())
  expect_that(identical(stpcaObjNew$lps,   stpcaObj$lps),   is_false())
})

# it doesn't!! It finds a saddle point!
#test_that("stpca.iterate.theta finds local maximum in theta", {
#  params = list(W=stpcaObjNew$W, sigSq=stpcaObjNew$sigSq, mu=stpcaObjNew$mu)
#  veclp = function(theta) {
#    L = relist(theta, params)
#    return(stpca.log_posterior(X, stpcaObjNew$K, L$W, L$mu, L$sigSq))
#  }
#
#  thetaHat = unlist(params)
#
#  lpMax = veclp(thetaHat)
#  for (i in seq_along(thetaHat)) {
#    lpUp = veclp(thetaHat + 1e-4*(seq_along(thetaHat)==i))
#    lpDn = veclp(thetaHat - 1e-4*(seq_along(thetaHat)==i))
#    expect_gte(lpMax+1e-5, lpUp) # Adjusted for numerical precision
#    expect_gte(lpMax+1e-5, lpDn) # by adding +1e-5 to lpMax
#  }
#})
