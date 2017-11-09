for (cfi in 1:5) {
  if (cfi==1) {
    cov.fun.str = "cov.SE"
    cov.fun   = cov.SE
    cov.fun.d = cov.SE.d
    nHPs      = 2
  } else if (cfi==2) {
    cov.fun.str = "cov.independent"
    cov.fun   = cov.independent
    cov.fun.d = cov.independent.d
    nHPs      = 1
  } else if (cfi==3) {
    cov.fun.str = "cov.RQ"
    cov.fun   = cov.RQ
    cov.fun.d = cov.RQ.d
    nHPs      = 3
  } else if (cfi==4) {
    cov.fun.str = "cov.noisy.SE"
    cov.fun   = cov.noisy.SE
    cov.fun.d = cov.noisy.SE.d
    nHPs      = 3
  } else {
    cov.fun.str = "cov.noisy.RQ"
    cov.fun   = cov.noisy.RQ
    cov.fun.d = cov.noisy.RQ.d
    nHPs      = 4
  }

  context(cov.fun.str)
  for (replicate in 1:5) {
    beta = rnorm(nHPs)
    K    = cov.fun(  locations, beta=beta)
    dK   = cov.fun.d(locations, beta=beta)
    D    = distanceMatrix(locations)

    test_that(paste(cov.fun.str,"returns valid covariance matrix"), {
      expect_is(K, "Matrix")
      expect_that(all(is.finite(K)), is_true())
      expect_that(isSymmetric(K), is_true())
    })

    test_that(paste(cov.fun.str,"gradient matches numerical gradient"), {
      jac.num = jacobian(function(beta_) {
        c(as.matrix(cov.fun(locations, beta=beta_)))
      }, x=beta)
      jac.analytic = vapply(dK, function(dKi) c(as.matrix(dKi)), numeric(d*d))
      expect_equal(jac.analytic, jac.num)
    })

    test_that(paste(cov.fun.str,"can accept precomuted distance matrix"), {
      K2 = cov.fun(locations, beta=beta, D=D)
      expect_equal(K, K2)
    })

    test_that(paste(cov.fun.str,"can accept precomputed sparse distance matrix, giving sparse result"), {
      Ds = distanceMatrix(locations, max.dist=1.2)
      Ks = cov.fun(locations, beta=beta, D=Ds)
      supported = which(D<1.2)
      expect_equal(Ks[supported], K[supported])
      expect_equal(Ks[-supported], rep(0, d*d-length(supported)))
    })

    test_that(paste(cov.fun.str,"derivatives can accept precomuted distance matrix"), {
      dK2 = cov.fun.d(locations, beta=beta, D=D)
      expect_equal(dK, dK2)
    })

    test_that(paste(cov.fun.str,"derivatives can accept precomputed sparse distance matrix, giving sparse result"), {
      Ds  = distanceMatrix(locations, max.dist=1.2)
      dKs = cov.fun.d(locations, beta=beta, D=Ds)
      supported = which(D<1.2)
      for (i in 1:length(dKs)) {
        expect_equal(dKs[[i]][supported], dK[[i]][supported])
        expect_equal(dKs[[i]][-supported], rep(0, d*d-length(supported)))
      }
    })
  }
}
