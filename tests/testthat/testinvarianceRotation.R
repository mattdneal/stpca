test_that("Rotating W before E-step gives same expectations as rotating W and E simultaniously", {
    Xc = sweep(stpca$X, 2, stpca$muHat)

    # Build W and E, rotate both after E-step
    W  = matrix(rnorm(stpca$d*stpca$k), nrow=stpca$d, ncol=stpca$k)
    E1 = EM.E(Xc, W, stpca$sigSqHat)

    R = svd(W)$v %*% diag(sign(W[1,]))
    Wrot   = W %*% R
    Vmean1 = E1$Vmean %*% R
    Vvar1  = lapply(E1$Vvar, function(vv) t(R) %*% vv %*% R)

    # Build rotated W, E step gives auto-rotated E
    E2 = EM.E(Xc, Wrot, stpca$sigSqHat)
    Vmean2 = E2$Vmean
    Vvar2  = E2$Vvar

    # Check these two methods don't produce different results
    expect_equal(Vmean1, Vmean2)
    expect_equal(Vvar1,  Vvar2)
})
