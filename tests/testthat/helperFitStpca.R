set.seed(1)

## Generate data
library(functional)
library(Matrix)
library(numDeriv)

n     = 15
k     = 4
dim   = c(6, 6)
d     = prod(dim)
beta  = log(c(2, 0.2, 1e-4))
k_se  = functional::Curry(cov.noisy.SE, beta=beta)
sigSq = 1.8

dataset = synthesize_data_kern(n, k, dim, kern=k_se, noisesd=sqrt(sigSq))
locs = dataset$grid
X = as.matrix(dataset$X)

beta0 = c(cov.SE.beta0(X, locs, k), log(1e-3))

stpca <- StpcaModel$new(X, k, beta0, locs, cov.noisy.SE, cov.noisy.SE.d, nIter=50)

#beta0 = log(c(1, 0.6, 0.6, 1e-3))
#stpca <- StpcaModel$new(X, k, beta0, locs, cov.noisy.MR, cov.noisy.MR.d, nIter=50)

constr = list(
  ineqA = rbind(diag(3), -diag(3)),
  ineqB = c(-log(c(0.01, 0.01, 1e-5)),
            log(c(100*crossprod(X[1,]),
                             5, 10)))
)
stpcaTuned = stpca$copy()$tune_beta(constraints=constr)
