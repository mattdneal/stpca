#' Performs StPCA
#'
#' @param X Data
#' @param k Latent dimensionality
#' @param locations the coordinates of each dimenion in X
#' @param covar.fn covariance function to generate K_beta
#' @param covar.fn.d gradient of the covariance function with respect to hyperparameters
#' @param beta0 initial hyperparameters
#' @param trace amount of reporting. 0=none, 1=low, 2=high
#' @param report_iter Number of iterations between reports
#' @param max.dist Maximum distance between features to consider
#' @param maxit.inner number of inner iterations
#' @param maxit.outer number of outer iterations
#' @return An \code{stpca} object.
#' @export
#' @include util.R
#' @include statistical-quantities.R
#' @include stpca.init.R
#' @include stpca.iterate.R
#' @import gsl
#' @import fields
#' @examples
#' library(fields)
#' data(ozone2)
#'
#' # Missing data: Replace missing values by column means
#' X = ozone2$y
#' for (col in 1:ncol(X)) {
#'   ind = is.na(X[,col])
#'   X[ind,col] = mean(X[,col], na.rm=TRUE)
#' }
#' X = X/sd(X) # Scale for numerical reasons
#'
#' locations = ozone2$lon.lat
#' locations = apply(locations, 2, function(col) (col-min(col))/(max(col)-min(col)))
#'
#' model.stpca = stpca(X, 3, locations, cov.SE, cov.SE.d, beta0=log(c(1, 0.5)),
#'                   maxit.inner=20, maxit.outer=3, trace=0)
stpca <- function(X, k, locations, covar.fn, covar.fn.d=NULL, beta0=c(),
                  trace=0, report_iter=10, max.dist=Inf,
                  maxit.inner=20, maxit.outer=5) {
  stpcaObj = stpca.init(X, k, locations, covar.fn,
                        covar.fn.d, beta0, trace, max.dist)
  stpcaObj = stpca.iterate(stpcaObj, trace, report_iter,
                           max.dist, maxit.inner, maxit.outer)
  return(stpcaObj)
}