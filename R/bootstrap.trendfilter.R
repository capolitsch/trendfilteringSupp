bootstrap.trendfilter <- function(x, y, lambda.opt, sigma = NULL, B = 1000, x.eval.grid = x, k = 2, 
                                  bootstrap.method = "nonparametric", alpha = 0.05, return.full.ensemble = F,
                                  max_iter = 250, obj_tol = 1e-06, mc.cores = 1){
  if ( is.null(sigma) ){
    data <- data.frame(x=x,y=y,wts=1)
  }else{
    data <- data.frame(x=x,y=y,wts=1/sigma^2)
  }
  
  if ( bootstrap.method == "nonparametric" ){
    if ( mc.cores == 1 ){
      tf.boot.ensemble <- matrix(unlist(replicate(B,tf.estimator(nonparametric.resampler(data), lambda.opt, k,
                                                               x.eval.grid, max_iter = max_iter, obj_tol = obj_tol))),
                               ncol = B)
    }else{
      par.func <- function(b){
        boot.tf.estimate <- tf.estimator(nonparametric.resampler(data), lambda.opt, k, x.eval.grid, 
                                         max_iter = max_iter, obj_tol = obj_tol)
        return(boot.tf.estimate)
      }
      tf.boot.ensemble <- matrix(unlist(parallel::mclapply(1:B, par.func, mc.cores = mc.cores)), ncol = B)
    }
  }
  
  if ( bootstrap.method == "parametric" ){
    tf.estimate <- tf.estimator(data, lambda.opt, k, x.eval.grid)
    data$tf.estimate <- tf.estimate
    if ( mc.cores == 1 ){
      tf.boot.ensemble <- matrix(unlist(replicate(B,tf.estimator(parametric.sampler(data), lambda.opt, k,
                                                                 x.eval.grid, max_iter = max_iter, obj_tol = obj_tol))),
                                 ncol = B)
    }else{
      par.func <- function(b){
        boot.tf.estimate <- tf.estimator(parametric.sampler(data), lambda.opt, k, x.eval.grid, 
                                         max_iter = max_iter, obj_tol = obj_tol)
        return(boot.tf.estimate)
      }
      tf.boot.ensemble <- matrix(unlist(parallel::mclapply(1:B, par.func, mc.cores = mc.cores)), ncol = B)
    }
  }
  
  if ( bootstrap.method == "wild" ){
    data$tf.estimate <- tf.estimator(data, lambda.opt, k, x.eval.grid)
    data$tf.residuals <- data$y - data$tf.estimate
    
    if ( mc.cores == 1 ){
      tf.boot.ensemble <- matrix(unlist(replicate(B,tf.estimator(wild.sampler(data), lambda.opt, k,
                                                               x.eval.grid, max_iter = max_iter, obj_tol = obj_tol))),
                               ncol = B)
    }else{
      par.func <- function(b){
        boot.tf.estimate <- tf.estimator(wild.sampler(data), lambda.opt, k, x.eval.grid, 
                                         max_iter = max_iter, obj_tol = obj_tol)
        return(boot.tf.estimate)
      }
      tf.boot.ensemble <- matrix(unlist(parallel::mclapply(1:B, par.func, mc.cores = mc.cores)), ncol = B)
    }
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
