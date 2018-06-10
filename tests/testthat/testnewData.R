context("Latent representations & reconstructions")

test_that("New data can be encoded to latent space", {
  Xnew <- stpcaUp$simulate(10, Wknown=TRUE)$X
  Vnew <- stpcaUp$encode(Xnew)$Vmean

  expect_equal(nrow(Xnew), nrow(Vnew))
  expect_true(all(is.finite(Vnew)))
})

test_that("Latent representation can be decoded into example", {
  Vnew <- matrix(rnorm(10*k), ncol=k)
  Xnew <- stpcaUp$decode(Vnew)

  expect_equal(nrow(Xnew), nrow(Vnew))
  expect_true(all(c(is.finite(Xnew))))
})

test_that("Larger magnitude latent repsresentations decode into larger samples", {
  Vnew     <- matrix(rnorm(10*k), ncol=k)
  VnewBig  <- Vnew * 100

  Xnew    <- sweep(stpcaUp$decode(Vnew),    2, stpcaUp$muHat)
  XnewBig <- sweep(stpcaUp$decode(VnewBig), 2, stpcaUp$muHat)

  expect_equal(c(XnewBig/Xnew), rep(100, length(Xnew)))
})

test_that("Going from latent to observed space and back to latent recovers analytic mean (shrunk toards 0)", {
  Vnew  <- matrix(rnorm(10*k), ncol=k)
  Xrec  <- stpcaUp$decode(Vnew)
  Vnew2 <- stpcaUp$encode(Xrec)$Vmean

  M <- crossprod(stpcaUp$WHat) + stpcaUp$sigSqHat*diag(k)
  shrinkFactor <- diag(k) - stpcaUp$sigSqHat*solve(M)
  expect_equal(Vnew2, Vnew%*%shrinkFactor)
})
