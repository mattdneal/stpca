context("Setters")

test_that("Can set new X", {
  ss <- stpcaUp$copy()

  Xnew <- (X*5 + 10)[1:10,]

  ss$set_X(Xnew)

  expect_equal(ss$X, Xnew)

  expect_error(ss$set_X(X[,1:10]))
})

test_that("Can set new locations", {
  ss <- stpcaUp$copy()
  Kold <- ss$K
  locsOld <- ss$locs

  # Shifting locs by 1 should not change K (b/c SE)
  ss$set_locs(ss$locs + 1)
  expect_equal(ss$K, Kold)
  expect_false(isTRUE(all.equal(ss$locs, locsOld)))

  # K should change here
  ss$set_locs(matrix(rnorm(ss$d*3), ncol=3))
  expect_equal(dim(Kold), dim(ss$K))
  expect_false(isTRUE(all.equal(Kold, ss$K)))

  expect_error(ss$set_locs(1))
})
