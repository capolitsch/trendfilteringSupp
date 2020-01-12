#' @keywords internal

parametric.sampler <- function(data){
  boot.sample <- data$tf.estimate + rnorm(nrow(data), sd = 1 / sqrt(data$wts))
  return(data.frame(x=data$x,y=boot.sample,wts=data$wts))
}
