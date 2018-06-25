context("Rotational invariance")

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

test_that("Log posterior is invariant to W -> WR", {
    K   <- stpcaUp$K
    W   <- stpcaUp$WHat
    sigSq <- stpcaUp$sigSqHat
    R1  <- svd(matrix(rnorm(k*k), ncol=k))$u
    R2  <- svd(matrix(rnorm(k*k), ncol=k))$v
    lp1 <- log_likelihood(Xc, W,      0, sigSq) + log_prior(K, W)
    lp2 <- log_likelihood(Xc, W%*%R1, 0, sigSq) + log_prior(K, W%*%R1)
    lp3 <- log_likelihood(Xc, W%*%R2, 0, sigSq) + log_prior(K, W%*%R2)

    expect_equal(lp1, lp2)
    expect_equal(lp1, lp3)
})

test_that("Complete log posterior is invariant to W -> WR, V -> VR", {
    K    <- stpcaUp$K
    W    <- stpcaUp$WHat
    V    <- stpcaUp$Vmean
    Vvar <- stpcaUp$Vvar
    sigSq <- stpcaUp$sigSqHat

    R1 <- svd(matrix(rnorm(k*k), ncol=k))$u
    R2 <- svd(matrix(rnorm(k*k), ncol=k))$v

    W1 <- W %*% R1
    W2 <- W %*% R2

    V1 <- V %*% R1
    V2 <- V %*% R2

    Vvar1 <- lapply(Vvar, function(vvar) {
      t(R1) %*% vvar %*% R1
    })

    Vvar2 <- lapply(Vvar, function(vvar) {
      t(R2) %*% vvar %*% R2
    })

    clp1 <- complete_log_posterior(Xc, V , Vvar , W , 0, sigSq, K)
    clp2 <- complete_log_posterior(Xc, V1, Vvar1, W1, 0, sigSq, K)
    clp3 <- complete_log_posterior(Xc, V2, Vvar2, W2, 0, sigSq, K)

    expect_equal(clp1, clp2)
    expect_equal(clp1, clp3)
})

## Checking if the sparse versions are invariant to R. They are NOT!!
## Analytically they DO change under W -> WR !
# test_that("Log posterior is invariant to W -> WR (sparse)", {
#     K   <- stpcaUp$K
#     W   <- stpcaUp$WHat
#     sigSq <- stpcaUp$sigSqHat
#     R1  <- svd(matrix(rnorm(k*k), ncol=k))$u
#     R2  <- svd(matrix(rnorm(k*k), ncol=k))$v
# 
#     lp1 <- log_likelihood(Xc, W,      0, sigSq) +
#       log_sparse_prior(W, K,      0.01)
# 
#     lp2 <- log_likelihood(Xc, W%*%R1, 0, sigSq) +
#       log_sparse_prior(W%*%R1, K, 0.01)
# 
#     lp3 <- log_likelihood(Xc, W%*%R2, 0, sigSq) +
#       log_sparse_prior(W%*%R2, K, 0.01)
# 
#     expect_equal(lp1, lp2)
#     expect_equal(lp1, lp3)
# })
# 
# test_that("Complete log posterior is invariant to W -> WR, V -> VR (sparse)", {
#     K    <- stpcaUp$K
#     W    <- stpcaUp$WHat
#     V    <- stpcaUp$Vmean
#     Vvar <- stpcaUp$Vvar
#     sigSq <- stpcaUp$sigSqHat
# 
#     R1 <- svd(matrix(rnorm(k*k), ncol=k))$u
#     R2 <- svd(matrix(rnorm(k*k), ncol=k))$v
# 
#     W1 <- W %*% R1
#     W2 <- W %*% R2
# 
#     V1 <- V %*% R1
#     V2 <- V %*% R2
# 
#     Vvar1 <- lapply(Vvar, function(vvar) {
#       t(R1) %*% vvar %*% R1
#     })
# 
#     Vvar2 <- lapply(Vvar, function(vvar) {
#       t(R2) %*% vvar %*% R2
#     })
# 
#     clp1 <- complete_log_posterior(Xc, V , Vvar , W , 0, sigSq, K, sparse=TRUE, b=0.01)
#     clp2 <- complete_log_posterior(Xc, V1, Vvar1, W1, 0, sigSq, K, sparse=TRUE, b=0.01)
#     clp3 <- complete_log_posterior(Xc, V2, Vvar2, W2, 0, sigSq, K, sparse=TRUE, b=0.01)
# 
#     expect_equal(clp1, clp2)
#     expect_equal(clp1, clp3)
# })
# 


