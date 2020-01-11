SURE.trendfilter <- function(x = NULL, y, sigma = NULL, lambda = NULL, k = 2, max_iter = 250, obj_tol = 1e-06){
  if ( is.null(x) ){
    x <- rep(1, length(y))
  }
  if ( is.null(sigma) ){
    stop("sigma must be provided in order to compute SURE.")
  }
  if ( length(sigma) == 1 ){
    sigma <- rep(sigma, length(y))
  }
  if ( is.null(lambda) ){
    stop("lambda must be specified.")
  }
  wts <- 1/sigma^2
  out <- glmgen::trendfilter(x = x, y = y, weights = wts, k = k, lambda = lambda,
                             control = glmgen::trendfilter.control.list(max_iter = max_iter, obj_tol = obj_tol))
  if ( length(lambda) == 1 ){
    SURE.loss <- mean( (out$beta - y)^2 ) + (2 * mean(1/wts) / length(x)) * out$df
  }
  if ( length(lambda) > 1 ){
    SURE.loss <- colMeans( (out$beta - y)^2 ) + (2 * mean(1/wts) / length(x)) * out$df
  }
  return(list(lambda = lambda, SURE.loss = as.numeric(SURE.loss)))
}
