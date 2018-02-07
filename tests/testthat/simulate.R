test_that("simulating new samples works ok", {
  Xsim <- stpcaUp$simulate(5, Wknown=FALSE)
  expect_equal(dim(Xsim), c(n, d))
  expect_true(all(is.finite(Xsim)))

  Xsim <- stpcaUp$simulate(5, Wknown=TRUE)
  expect_equal(dim(Xsim), c(n, d))
  expect_true(all(is.finite(Xsim)))
})
