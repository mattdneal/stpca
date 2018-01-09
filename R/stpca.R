#' Structured PCA Model
#'
#' @import Matrix
#' @import maxLik
#' @export
StpcaModel <- setRefClass("StpcaModel",
  fields = list(
    X        = "matrix",
    n        = "integer",
    k        = "integer",
    d        = "integer",
    beta     = "numeric",
    K        = "Matrix",
    KD       = "list",
    locs     = "matrix",
    covFn    = "function",
    covFnD   = "function",
    muHat    = "numeric",
    WHat     = "matrix",
    sigSqHat = "numeric",
    Vmean    = "matrix",
    Vvar     = "list",
    logEvidences  = "numeric",
    logPosteriors = "numeric",
    logEvidenceD  = "numeric"
  ),
  methods = list(
    initialize = function(X=matrix(nrow=0, ncol=0), k=1, beta0=numeric(0), locs=matrix(),
                          covFn=function() NULL, covFnD=function() NULL, ...) {
      X      <<- X
      n      <<- nrow(X)
      k      <<- as.integer(k)
      d      <<- ncol(X)
      locs   <<- locs
      covFn  <<- covFn
      covFnD <<- covFnD

      if (nrow(X)>0) {
        set_beta(beta0)
        thetaInit <- initialize_from_ppca(X, k)
        muHat    <<- thetaInit$mu
        sigSqHat <<- thetaInit$sigSq
        WHat     <<- thetaInit$W
      }

      callSuper(...)
    },
    fit = function(nIter=50) {
      vals <- fit_stpca(X, WHat, muHat, sigSqHat, K, nIter)

      WHat          <<- vals$WHat
      sigSqHat      <<- vals$sigSqHat
      Vmean         <<- vals$Vmean
      Vvar          <<- vals$Vvar
      logEvidences  <<- c(logEvidences,  vals$logEvidence)
      logPosteriors <<- c(logPosteriors, vals$logPosteriors)

      return(.self)
    },
    set_beta = function(beta, nIter=50) {
      'Documentations for the method goes here'
      beta <<- beta
      K    <<- covFn(locs, beta=beta)
      fit()
      return(.self)
    },
    compute_gradient = function() {
      KD <<- covFnD(locs, beta=beta)
      logEvidenceD <<- log_evidence_d(X, K, WHat, muHat, sigSqHat, beta, KD)
      return(.self)
    },
    tune_beta = function(nIter) {
      stop("Implement me!")
      return(.self)
    }
  )
)
