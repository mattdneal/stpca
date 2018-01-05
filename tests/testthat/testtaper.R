context("Covariance function tapering")
test_that("Tapering works sensibly", {
  K.SE   = cov.SE(locs, beta=c(2, 0.5))
  K.SE.t = cov.taper(cov.SE)(locs, beta=c(2, 0.5))
  expect_equal(diag(K.SE), diag(K.SE.t))

  upper1 = K.SE[upper.tri(K.SE)]
  upper2 = K.SE.t[which(upper.tri(K.SE.t))]

  expect_that(all(upper1 > upper2), is_true())

  expect_is(K.SE.t, "sparseMatrix")
})

test_that("Tapered gradients match numerical gradients", {
  beta = c(0.5, 1.2)

  jac.num = jacobian(function(beta_) {
    c(as.matrix(cov.taper(cov.SE)(locs, beta=beta_)))
  }, x=beta)

  KD = cov.taper.d(cov.SE.d)(locs, beta=beta)
  jac.analytic = vapply(KD, function(KDi) c(as.matrix(KDi)), numeric(d*d))
  expect_equivalent(jac.analytic, jac.num)
})

test_that("Tapered covariance functions can accept (and ignore) D", {
  require(fields)
  D = rdist(locs)
  expect_error(
    cov.taper(cov.SE)(locs, beta=rnorm(2), D=rdist(D)),
    NA # Assert that *no* error is expected
  )
  expect_error(
    cov.taper.d(cov.SE.d)(locs, beta=rnorm(2), D=rdist(D)),
    NA
  )
})
