set.seed(1)

## Generate data
library(functional)
library(Matrix)
library(numDeriv)

n     = 15
k     = 4
dim   = c(5, 5)
d     = prod(dim)
beta  = log(c(2, 0.2, 1e-4))
k_se  = functional::Curry(cov.noisy.SE, beta=beta)
sigSq = 1.8

dataset = synthesize_data_kern(n, k, dim, kern=k_se, noisesd=sqrt(sigSq))
locs = dataset$grid
X = as.matrix(dataset$X)

beta0 = c(cov.SE.beta0(X, locs, k), log(1e0))

stpca <- StpcaModel$new(X, k, beta0, locs, cov.noisy.SE, cov.noisy.SE.d, maxit=100)

#beta0 = log(c(1, 0.6, 0.6, 1e-3))
#stpca <- StpcaModel$new(X, k, beta0, locs, cov.noisy.MR, cov.noisy.MR.d, maxit=50)

constr = list(
  ineqA = rbind(diag(3), -diag(3)),
  ineqB = c(-log(c(0.01, 0.01, 1e-5)),
            log(c(100*crossprod(X[1,]),
                             5, 2)))
)

stpcaUpBeta  = stpca$copy()$update_beta(constraints=constr)
stpcaUpTheta = stpcaUpBeta$copy()$update_theta(100, bftol=1e-6)

# Set beta to value near max so fewer updates can be done to converge.
goodBeta = c(-0.965644, -1.949905, -6.886124)
stpcaUp = stpca$
            copy()$
            set_beta(goodBeta)$
            update_theta()$
            update(5, 100, constraints=constr)
