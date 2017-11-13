set.seed(1)

## Generate data
library(functional)
library(Matrix)
library(numDeriv)

n     = 15
k     = 4
dim   = c(5, 5)
d     = prod(dim)
beta  = log(c(2, 0.2))
k_se  = functional::Curry(cov.SE, beta=beta)
sigSq = 1.8

dataset = synthesize_data_kern(n, k, dim, kern=k_se, noisesd=sqrt(sigSq))
locations = dataset$grid
X = as.matrix(dataset$X)

beta0 = c(cov.SE.beta0(X, locations, k), log(1e-3))

stpcaObj = stpca(X, k, locations, cov.noisy.SE, cov.noisy.SE.d, beta0, trace=0, maxit.inner=0, maxit.outer=0)

H = stpca.H(X, stpcaObj$W, stpcaObj$mu, stpcaObj$sigSq, stpcaObj$K)

maxit.inner=20
maxit.outer=1

stpcaObj.it = stpca.iterate(stpcaObj, maxit.inner=maxit.inner,
                                      maxit.outer=maxit.outer)

stpcaObj.it.b = stpca.iterate.beta(stpcaObj.it)
stpcaObj.it.t = stpca.iterate.theta(stpcaObj.it.b, maxit.inner=50)
