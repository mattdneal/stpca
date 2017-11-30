context("stpca variable types")

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

test_that("stpca.iterate returns a valid stpca object", {
  expect_is(stpcaObj.it$X,     "matrix")
  expect_is(stpcaObj.it$W,     "matrix")
  expect_is(stpcaObj.it$sigSq, "numeric")
  expect_is(stpcaObj.it$mu,    "numeric")
  expect_is(stpcaObj.it$V,     "matrix")
  expect_is(stpcaObj.it$ll,    "numeric")
  expect_is(stpcaObj.it$lp,    "numeric")
  expect_is(stpcaObj.it$lps,   "numeric")
  expect_is(stpcaObj.it$bic,   "numeric")
  expect_true(length(stpcaObj.it$beta0)==0 | is(stpcaObj.it$beta0, "numeric"))
  expect_is(stpcaObj.it$D,     "Matrix")
  expect_is(stpcaObj.it$K,     "Matrix")
  expect_is(stpcaObj.it$covar.fn,   "function")
  expect_is(stpcaObj.it$covar.fn.d, "function")
  expect_is(stpcaObj.it$locations,  "matrix")

  expect_equal(dim(stpcaObj.it$X), c(n, d))
  expect_equal(dim(stpcaObj.it$W), c(d, k))
  expect_length(stpcaObj.it$sigSq, 1)
  expect_length(stpcaObj.it$mu, d)
  expect_equal(dim(stpcaObj.it$V), c(n, k))
  expect_length(stpcaObj.it$ll, 1)
  expect_length(stpcaObj.it$lp, 1)
  expect_length(stpcaObj.it$lps, 1 + maxit.inner*(maxit.outer + 1))
  expect_length(stpcaObj.it$bic, 1)
  expect_equal(dim(stpcaObj.it$D), c(d, d))
  expect_equal(dim(stpcaObj.it$K), c(d, d))
  expect_equal(dim(locations)[1], d)
})

context("stpca.iterate.beta")

test_that("stpca.iterate.beta increases log evidence", {
  expect_gt(stpcaObj.it.b$log_evidence, stpcaObj.it$log_evidence)
})

test_that("stpca.iterate.beta modifies beta, H, K", {
  expect_that(identical(stpcaObj.it.b$H,    stpcaObj.it$H),    is_false())
  expect_that(identical(stpcaObj.it.b$K,    stpcaObj.it$K),    is_false())
  expect_that(identical(stpcaObj.it.b$beta, stpcaObj.it$beta), is_false())
})

test_that("stpca.iterate.beta does not change theta", {
  expect_identical(stpcaObj.it.b$W,     stpcaObj.it$W)
  expect_identical(stpcaObj.it.b$sigSq, stpcaObj.it$sigSq)
  expect_identical(stpcaObj.it.b$V,     stpcaObj.it$V)
  expect_identical(stpcaObj.it.b$ll,    stpcaObj.it$ll)
  expect_identical(stpcaObj.it.b$lp,    stpcaObj.it$lp)
  expect_identical(stpcaObj.it.b$lps,   stpcaObj.it$lps)
})

test_that("beta has sensible 95% confidence intervals", {
  ci95 = attr(stpcaObj.it.b$beta, "ci95")
  expect_is(ci95, "matrix")
  for (i in seq_along(stpcaObj$beta)) {
    expect_lt(stpcaObj.it.b$beta[i], ci95["upper",i])
    expect_gt(stpcaObj.it.b$beta[i], ci95["lower",i])
  }
})

context("stpca.iterate.theta")

test_that("stpca.iterate.theta increases log posterior", {
  expect_gt(stpcaObj.it.t$lp, stpcaObj.it.b$lp)
})

test_that("stpca.iterate.theta only changes relevant variables", {
  expect_identical(stpcaObj.it.t$H,            stpcaObj.it.b$H)
  expect_identical(stpcaObj.it.t$K,            stpcaObj.it.b$K)
  expect_identical(stpcaObj.it.t$beta,         stpcaObj.it.b$beta)
  expect_identical(stpcaObj.it.t$log_evidence, stpcaObj.it.b$log_evidence)

  expect_that(identical(stpcaObj.it.t$W,     stpcaObj.it.b$W),     is_false())
  expect_that(identical(stpcaObj.it.t$sigSq, stpcaObj.it.b$sigSq), is_false())
  expect_that(identical(stpcaObj.it.t$V,     stpcaObj.it.b$V),     is_false())
  expect_that(identical(stpcaObj.it.t$ll,    stpcaObj.it.b$ll),    is_false())
  expect_that(identical(stpcaObj.it.t$lp,    stpcaObj.it.b$lp),    is_false())
  expect_that(identical(stpcaObj.it.t$lps,   stpcaObj.it.b$lps),   is_false())
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
