library(devtools)
install_github("statsmaths/glmgen", subdir="R_pkg/glmgen")

Cp_trendfilter <- function(x = NULL, y, wts = NULL, deg = 1, lambda = NULL){
  if ( is.null(x) ){
    x <- rep(1, length(y))
  }
  if ( is.null(wts) ){
    stop("Weights must be provided in order to compute Mallows' Cp.")
  }
  if ( is.null(lambda) ){
    stop("lambda must be specified.")
  }
  out <- glmgen::trendfilter(x = x, y = y, weights = wts, k = deg, lambda = lambda)
  if ( length(lambda) == 1 ){
    Cp.loss <- mean( (out$beta - y)^2 ) + (2 * mean(1/wts) / length(x)) * out$df
  }
  if ( length(lambda) > 1 ){
    Cp.loss <- colMeans( (out$beta - y)^2 ) + (2 * mean(1/wts) / length(x)) * out$df
  }
  return(list(lambda = lambda, Cp.loss = as.numeric(Cp.loss)))
}

# Example
lambda.grid <- exp(seq(-16,2,length=50))
