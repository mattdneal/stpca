for (cfi in 1:6) {
  if (cfi==1) {
    cov.fun.str = "cov.SE"
    cov.fun   = cov.SE
    cov.fun.d = cov.SE.d
    nHPs      = 2
    acceptD = TRUE
  } else if (cfi==2) {
    cov.fun.str = "cov.independent"
    cov.fun   = cov.independent
    cov.fun.d = cov.independent.d
    nHPs      = 1
    acceptD = TRUE
  } else if (cfi==3) {
    cov.fun.str = "cov.RQ"
    cov.fun   = cov.RQ
    cov.fun.d = cov.RQ.d
    nHPs      = 3
    acceptD = TRUE
  } else if (cfi==4) {
    cov.fun.str = "cov.noisy.SE"
    cov.fun   = cov.noisy.SE
    cov.fun.d = cov.noisy.SE.d
    nHPs      = 3
    acceptD = TRUE
  } else if (cfi==5) {
    cov.fun.str = "cov.noisy.RQ"
    cov.fun   = cov.noisy.RQ
    cov.fun.d = cov.noisy.RQ.d
    nHPs      = 4
    acceptD = TRUE
  } else if (cfi==6) {
    cov.fun.str = "cov.MR"
    cov.fun   = cov.MR
    cov.fun.d = cov.MR.d
    nHPs      = ncol(locs)+1
    acceptD = FALSE
  } else {
    cov.fun.str = "cov.noisy.MR"
    cov.fun   = cov.noisy.MR
    cov.fun.d = cov.noisy.MR.d
    nHPs      = ncol(locs)+2
    acceptD = FALSE
  }

  context(cov.fun.str)
  # 6 Replicates for different values of beta.
  for (replicate in 1:6) {
    beta = rnorm(nHPs)
    K    = cov.fun(  locs, beta=beta)
    KD   = cov.fun.d(locs, beta=beta)
    D    = distanceMatrix(locs)

    test_that(paste(cov.fun.str,"returns valid covariance matrix"), {
      expect_is(K, "Matrix")
      expect_that(all(is.finite(K)), is_true())
      expect_that(isSymmetric(K), is_true())
      expect_is(K, "symmetricMatrix")
    })

    test_that(paste(cov.fun.str,"gradient matches numerical gradient"), {
      jac.num = jacobian(function(beta_) {
        c(as.matrix(cov.fun(locs, beta=beta_)))
      }, x=beta)
      jac.analytic = vapply(KD, function(KDi) c(as.matrix(KDi)), numeric(d*d))
      expect_equivalent(jac.analytic, jac.num)
    })

    if (acceptD) {
      test_that(paste(cov.fun.str,"can accept precomputed distance matrix"), {
        K2 = cov.fun(locs, beta=beta, D=D)
        expect_equal(K, K2)
      })

      test_that(paste(cov.fun.str,"can accept precomputed sparse distance matrix, giving sparse result"), {
        Ds = distanceMatrix(locs, max.dist=1.2)
        Ks = cov.fun(locs, beta=beta, D=Ds)
        supported = which(D<1.2)
        expect_equal(Ks[supported], K[supported])
        expect_equal(Ks[-supported], rep(0, d*d-length(supported)))
      })

      test_that(paste(cov.fun.str,"derivatives can accept precomuted distance matrix"), {
        KD2 = cov.fun.d(locs, beta=beta, D=D)
        expect_equal(KD, KD2)
      })

      test_that(paste(cov.fun.str,"derivatives can accept precomputed sparse distance matrix, giving sparse result"), {
        Ds  = distanceMatrix(locs, max.dist=1.2)
        KDs = cov.fun.d(locs, beta=beta, D=Ds)
        supported = which(D<1.2)
        for (i in 1:length(KDs)) {
          expect_equal(KDs[[i]][supported], KD[[i]][supported])
          expect_equal(KDs[[i]][-supported], rep(0, d*d-length(supported)))
        }
      })
    }
  }
}
