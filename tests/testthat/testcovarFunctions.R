context("covariance functions")

library(Matrix)

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

test_that("cov.SE can accept precomputed distance matrix", {
  set.seed(1)
  beta = rnorm(2)
  K1 = cov.SE(locations, beta=beta)

  D = distanceMatrix(locations) # Nonsparse: explicit zeroes
  K2 = cov.SE(locations, beta=beta, D=D)

  expect_equal(K1, K2)
})

test_that("cov.independant works as expected", {
  set.seed(1)
  beta = rnorm(1)
  K = cov.independent(locations, beta=beta)
  expect_is(K, "Matrix")
  expect_that(all(is.finite(K)), is_true())
  expect_equivalent(K, Diagonal(d, exp(beta)))
})

test_that("cov.independant.d matches numerical gradient", {
  set.seed(1)
  beta = rnorm(1)
  dK = cov.independent.d(locations, beta=beta)

  jac.num = c(jacobian(function(beta_) {
    c(as.matrix(cov.independent(locations, beta=beta_)))
  }, x=beta))
  jac.analytic = c(as.matrix(dK))
  expect_equivalent(jac.analytic, jac.num)
})

test_that("cov.independant can accept precomputed distance matrix", {
  set.seed(1)
  beta = rnorm(1)
  K1 = cov.independent(locations, beta=beta)

  D = distanceMatrix(locations) # Nonsparse: explicit zeroes
  K2 = cov.independent(locations, beta=beta, D=D)

  expect_equal(K1, K2)

  D = distanceMatrix(locations, max.dist=0.5) # Sparse
  K3 = cov.independent(locations, beta=beta, D=D)

  expect_equal(K1, K3)
})
