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
    logEvidence   = "numeric",
    logPosteriors = "numeric",
    logEvidenceD  = "numeric",
    H        = "list",
    maxim    = "maxim"
  ),
  methods = list(
    initialize = function(X=matrix(nrow=0, ncol=0), k=1, beta0=numeric(0),
                          locs=matrix(), covFn=function() NULL,
                          covFnD=function() NULL, nIter=50, ...) {
      X      <<- X
      n      <<- nrow(X)
      k      <<- as.integer(k)
      d      <<- ncol(X)
      locs   <<- locs
      covFn  <<- covFn
      covFnD <<- covFnD

      # Initialize the 'maxim' slot with an empty maxim object.
      emptyMaxim <- list()
      class(emptyMaxim) <- "maxim"
      maxim  <<- emptyMaxim

      if (nrow(X)>0) {
        thetaInit <- initialize_from_ppca(X, k)
        muHat    <<- thetaInit$mu
        sigSqHat <<- thetaInit$sigSq
        WHat     <<- thetaInit$W
        set_beta(beta0, nIter)
      }

      callSuper(...)
    },
    fit = function(nIter=50) {
      vals <- NULL
      try(vals <- fit_stpca(X, WHat, muHat, sigSqHat, K, nIter))

      if (is.null(vals)) {
        Vmean         <<- matrix()
        Vvar          <<- list()
        logEvidence   <<- -Inf
        logPosteriors <<- -Inf
        H             <<- list()
      } else {
        WHat          <<- vals$WHat
        sigSqHat      <<- vals$sigSqHat
        Vmean         <<- vals$Vmean
        Vvar          <<- vals$Vvar
        logEvidence   <<- vals$logEvidence
        logPosteriors <<- vals$logPosteriors
        H             <<- vals$H
      }

      invisible(.self)
    },
    set_beta = function(betaNew, nIter=50) {
      'Documentations for the method goes here'
      if (!identical(betaNew, beta)) {
        beta <<- betaNew
        K    <<- Matrix()
        try(K <<- covFn(locs, beta=betaNew), silent=TRUE) # As dppMatrix?

        if (identical(K, Matrix()) || any(!is.finite(K@x))) {
          # Enter this branch if K is problematic. If K could not be built,
          # then it defaults to Matrix(). If K could be build but is bad,
          # it will contain non-finite values.
          Vmean         <<- matrix()
          Vvar          <<- list()
          logEvidence   <<- -Inf
          logPosteriors <<- -Inf
        } else {
          fit()
        }
      }
      invisible(.self)
    },
    compute_gradient = function() {
      KD           <<- covFnD(locs, beta=beta)
      logEvidenceD <<- log_evidence_d(X, K, WHat, muHat, sigSqHat, beta, KD, H)
      invisible(.self)
    },
    tune_beta = function(...) {
      logLik = function(beta_) {
        set_beta(beta_)$logEvidence
      }

      logLikGrad = function(beta_) {
        set_beta(beta_)$compute_gradient()$logEvidenceD
      }

      maxim <<- maxBFGS(logLik, logLikGrad, start=beta, ...)
      set_beta(maxim$estimate)

      invisible(.self)
    }
  )
)
