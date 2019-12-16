bootstrap.trendfilter <- function(x, y, lambda.opt, sigma = NULL, B = 1000, x.eval.grid = x, k = 2, 
                                  bootstrap.method = "nonparametric", alpha = 0.05, return.full.ensemble = F, max_iter = 250){
  if ( is.null(sigma) ){
    data <- data.frame(x=x,y=y,wts=1)
  }else{
    data <- data.frame(x=x,y=y,wts=1/sigma^2)
  }
  
  if ( bootstrap.method == "nonparametric" ){
    tf.boot.ensemble <- matrix(unlist(replicate(B,tf.estimator(nonparametric.resampler(data), lambda.opt, k, x.eval.grid,
                                                              control = trendfilter.control.list(max_iter = max_iter)))),
                               ncol = B)
  }
  
  if ( bootstrap.method == "parametric" ){
    tf.estimate <- tf.estimator(data, lambda.opt, k, x.eval.grid)
    data$tf.estimate <- tf.estimate
    tf.boot.ensemble <- matrix(unlist(replicate(B,tf.estimator(parametric.sampler(data), lambda.opt, k, x.eval.grid,
                                                              control = trendfilter.control.list(max_iter = max_iter)))),
                               ncol = B)
  }
  
  if ( bootstrap.method == "wild" ){
    tf.estimate <- tf.estimator(data, lambda.opt, k, x.eval.grid)
    tf.residuals <- data$y - tf.estimate
    data$tf.estimate <- tf.estimate
    data$tf.residuals <- tf.residuals
    tf.boot.ensemble <- matrix(unlist(replicate(B,tf.estimator(wild.sampler(data), lambda.opt, k, x.eval.grid,
                                                              control = trendfilter.control.list(max_iter = max_iter)))),
                               ncol = B)
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
