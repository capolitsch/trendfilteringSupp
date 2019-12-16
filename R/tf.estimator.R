tf.estimator <- function(data, lambda, k, x.eval.grid){
  tf.fit <- glmgen::trendfilter(data$x, data$y, data$wts, k = 2, lambda = lambda)
  tf.estimate <- as.numeric(suppressWarnings(glmgen:::predict.trendfilter(tf.fit, x.new = x.grid, lambda = lambda)))
  return(tf.estimate)
}