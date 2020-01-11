wild.sampler <- function(data){
  wild.boot.resids <- data$tf.residuals * sample(x = c((1+sqrt(5))/2, 1-sqrt(5)/2), size = nrow(data), replace = T, 
                                                 prob = c((1+sqrt(5))/(2*sqrt(5)),(sqrt(5)-1)/(2*sqrt(5))))
  wild.boot.sample <- data$tf.estimate + wild.boot.resids
  return(data.frame(x=data$x,y=wild.boot.sample,wts=data$wts))
}
