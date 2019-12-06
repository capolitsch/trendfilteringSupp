install.packages("devtools")
devtools::install_github("statsmaths/glmgen", subdir="R_pkg/glmgen")
install.packages("dplyr")

#####################################################################################
# Uncertainty quantification with trend filtering via various bootstrap methods
#####################################################################################
# bootstrap.trendfilter takes data of the form (x_1,y_1), ...,(x_n,y_n) or
# (x_1,y_1,sigma_1), ...,(x_n,y_n,sigma_n) and implements one of three bootstrap algorithms
# 
# The bootstrap.method should generally be chosen according to the following criteria:
#
# S1. The inputs are irregularly sampled --> bootstrap.method = "nonparametric"
# S2. The inputs are regularly sampled and the noise distribution is known --> bootstrap.method = "parametric"
# S3. The inputs are regularly sampled and the noise distribution is unknown --> bootstrap.method = "wild"
#
# lambda.opt is the optimally chosen trend filtering hyperparameter obtained by, e.g. minimizing SURE (see SURE.trendfilter),
# cross validation (see glmgen::cv.trendfilter), or chi-squared (see chi.squared.trendfilter)
#
# bootstrap.trendfilter outputs the 1-alpha pointwise variability bands obtained by computing the sample quantiles of
# the bootstrap ensemble
#
# if return.full.ensemble = TRUE, the full bootstrap ensemble will also be returned as a matrix

nonpar.resampler <- function(data){
  resampled.data <- dplyr::sample_n(data, size = nrow(data))
  return(resampled.data)
}

parametric.sampler <- function(data){
  boot.sample <- data$tf.estimate + rnorm(nrow(data), sd = 1 / sqrt(data$wts))
  return(boot.sample)
}

wild.sampler <- function(data){
  wild.boot.resids <- data$tf.residuals * sample(x = c((1+sqrt(5))/2, 1-sqrt(5)/2), size = nrow(data), replace = T, 
                                                 prob = c((1+sqrt(5))/(2*sqrt(5)),(sqrt(5)-1)/(2*sqrt(5))))
  wild.boot.sample <- data$tf.estimate + wild.boot.resids
  return(wild.boot.sample)
}

tf.estimator <- function(data, lambda, k, x.eval.grid){
  tf.fit <- glmgen::trendfilter(data$x, data$y, data$wts, k = 2, lambda = lambda)
  tf.estimate <- as.numeric(suppressWarnings(glmgen:::predict.trendfilter(tf.fit, x.new = x.grid, lambda = lambda)))
  return(tf.estimate)
}

bootstrap.trendfilter <- function(x, y, lambda.opt, sigma = NULL, B = 1000, x.eval.grid = x, k = 2, 
                                  bootstrap.method = "nonparametric", alpha = 0.05, return.full.ensemble = F){
  if ( is.null(sigma) ){
    data <- data.frame(x=x,y=y,wts=1/sigma^2)
  }else{
    data <- data.frame(x=x,y=y,wts=1)
  }
  
  if ( bootstrap.method == "nonparametric" ){
    tf.boot.ensemble <- matrix(unlist(replicate(B,tf.estimator(nonpar.resampler(data), lambda.opt, k, x.eval.grid))), ncol = B)
  }
  
  if ( bootstrap.method == "parametric" ){
    tf.estimate <- tf.estimator(data, lambda.opt, k, x.eval.grid)
    data$tf.estimate <- tf.estimate
    tf.boot.ensemble <- matrix(unlist(replicate(B,tf.estimator(parametric.sampler(data), lambda.opt, k, x.eval.grid))), ncol = B)
  }
  
  if ( bootstrap.method == "wild" ){
    tf.estimate <- tf.estimator(data, lambda.opt, k, x.eval.grid)
    tf.residuals <- data$y - tf.estimate
    data$tf.estimate <- tf.estimate
    data$tf.residuals <- tf.residuals
    tf.boot.ensemble <- matrix(unlist(replicate(B,tf.estimator(wild.sampler(data), lambda.opt, k, x.eval.grid))), ncol = B)
  }
  
  bootstrap.lower.perc.intervals <- apply(tf.boot.ensemble,1,quantile,probs = alpha/2)
  bootstrap.upper.perc.intervals <- apply(tf.boot.ensemble,1,quantile,probs = 1-alpha/2)
  
  if ( !return.full.ensemble ){
    return(list(bootstrap.lower.perc.intervals=bootstrap.lower.perc.intervals,
                bootstrap.upper.perc.intervals=bootstrap.upper.perc.intervals))
  }else{
    return(list(bootstrap.lower.perc.intervals=bootstrap.lower.perc.intervals,
                bootstrap.upper.perc.intervals=bootstrap.upper.perc.intervals,
                tf.boot.ensemble = tf.boot.ensemble))
  }
}
