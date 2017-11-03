context("covariance functions")

locations = expand.grid(seq(-1, 1, length=10), seq(-1, 1, length=10))
beta = log(c(1.5, 0.9))
d = nrow(locations)

test_that("cov.SE works as expected", {
  K = cov.SE(locations, beta=beta)
  expect_is(K, "Matrix")
  expect_that(all(is.finite(K)), is_true())
  expect_equal(diag(K), rep(exp(beta[1]), nrow(locations)))
  expect_lt(max(K[upper.tri(K, diag=FALSE)]), exp(beta[1]))
})

test_that("cov.SE.d matches numerical gradient", {
  dK = cov.SE.d(locations, beta=beta)

  # Derivative wrt variance at zero distance should always be exp(beta[1])
  expect_equal(diag(dK[[1]]), rep(exp(beta[1]), d))

  # Derivative wrt lengthscale at zero distance should always be 0
  expect_equal(diag(dK[[2]]), rep(0, d))

  jac.num = jacobian(function(beta_) {
    c(as.matrix(cov.SE(locations, beta=beta_)))
  }, x=beta)
  jac.analytic = vapply(dK, function(dKi) c(as.matrix(dKi)), numeric(d*d))
  expect_equal(jac.analytic, jac.num)

  beta=rnorm(2)
  dK = cov.SE.d(locations, beta=beta)
  jac.num = jacobian(function(beta_) {
    c(as.matrix(cov.SE(locations, beta=beta_)))
  }, x=beta)
  jac.analytic = vapply(dK, function(dKi) c(as.matrix(dKi)), numeric(d*d))
  expect_equal(jac.analytic, jac.num)
})
