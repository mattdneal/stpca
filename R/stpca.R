#' Structured PCA Model
#'
#' @export
StpcaModel <- setRefClass("StpcaModel",
  fields = list(
    X        = "matrix",
    n        = "integer",
    k        = "integer",
    d        = "integer",
    beta     = "numeric",
    K        = "Matrix",
    locs     = "matrix",
    covFn    = "function",
    covFnD   = "function",
    muHat    = "numeric",
    WHat     = "matrix",
    sigSqHat = "numeric",
    Vmean    = "matrix",
    Vvar     = "list",
    logEvidences  = "numeric",
    logPosteriors = "numeric"
  ),
  methods = list(
    initialize = function(X=matrix(), k=1, beta0=numeric(0), locs=matrix(),
                          covFn=NULL, covFnD=NULL, ...) {
      X      <<- X
      n      <<- nrow(X)
      k      <<- as.integer(k)
      d      <<- ncol(X)
      beta   <<- beta0
      locs   <<- locs
      covFn  <<- covFn
      covFnD <<- covFnD

      ## Initialize theta (W, mu, sigSq) from PPCA
      Xc <- scale(X, scale=FALSE) # Centered data: nxd
      covar.svd <- svd(Xc/sqrt(n), nu=0, nv=k)
      covar.eigval <- covar.svd$d^2

      muHat    <<- attr(Xc, "scaled:center")
      sigSqHat <<- sum(covar.eigval[-(1:k)])/(d-k)
      WHat     <<- (covar.svd$v %*%
                    diag(sqrt(covar.eigval[seq_len(k)] - sigSqHat)))

      set_beta(beta0)

      callSuper(...)
    },
    set_beta = function(beta, nIter=50) {
      'Documentations for the method goes here'
      beta <<- beta
      K    <<- covFn(locs, beta=beta)

      vals <- fit_stpca(X, WHat, muHat, sigSqHat, K, nIter)

      WHat          <<- vals$WHat
      sigSqHat      <<- vals$sigSqHat
      Vmean         <<- vals$Vmean
      Vvar          <<- vals$Vvar
      logEvidences  <<- c(logEvidences,  vals$logEvidence)
      logPosteriors <<- c(logPosteriors, vals$logPosteriors)

      return(.self)
    },
    tune_beta = function(nIter) {
      stop("Implement me!")
      return(.self)
    }
  )
)
