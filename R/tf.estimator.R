tf.estimator <- function(data, lambda, k, x.eval.grid, max_iter = 250){
  tf.fit <- glmgen::trendfilter(data$x, data$y, data$wts, k = 2, lambda = lambda, trendfilter.control.list(max_iter = max_iter))
  tf.estimate <- as.numeric(suppressWarnings(glmgen:::predict.trendfilter(tf.fit, x.new = x.eval.grid, lambda = lambda)))
  return(tf.estimate)
}
