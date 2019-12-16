parametric.sampler <- function(data){
  boot.sample <- data$tf.estimate + rnorm(nrow(data), sd = 1 / sqrt(data$wts))
  return(boot.sample)
}