tf.estimator <- function(data, lambda, k, x.eval.grid, max_iter = 250, obj_tol = 1e-05){
  tf.fit <- glmgen::trendfilter(data$x, data$y, data$wts, k = 2, lambda = lambda, 
                                control = glmgen::trendfilter.control.list(max_iter = max_iter, obj_tol = obj_tol))
  tf.estimate <- as.numeric(suppressWarnings(glmgen:::predict.trendfilter(tf.fit, x.new = x.eval.grid, lambda = lambda)))
  return(tf.estimate)
}
