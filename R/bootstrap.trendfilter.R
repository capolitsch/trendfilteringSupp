#' Optimize the trend filtering hyperparameter with respect to Stein's unbiased 
#' risk estimate
#'
#' @description \code{bootstrap.trendfilter} implements one of three bootstrap 
#' algorithms to obtain variability bands to accompany an optimized trend
#' filtering point estimate. The bootstrap method should generally be chosen 
#' according to the following criteria: \cr \cr
#' S1. The inputs are irregularly 
#' sampled –> \code{bootstrap.method = "nonparametric"} \cr \cr
#' S2. The inputs are regularly 
#' sampled and the noise distribution is known –> \code{bootstrap.method = "parametric"} \cr \cr
#' S3. The inputs are regularly sampled and the noise distribution is unknown –> 
#' \code{bootstrap.method = "wild"} \cr
#' @param x A vector of the observed inputs. If \code{NULL}, then we assume
#' unit spacing.
#' @param y A vector of the observed outputs.
#' @param sigma A vector of measurement standard errors for the observed outputs.
#' @param lambda.min The optimally-tuned trend filtering hyperparameter, e.g. by 
#' minimizing SURE (see SURE.trendfilter) or cross validation.
#' @param k The degree of the trend filtering estimator. Defaults to \code{k=2}
#' (quadratic trend filtering).
#' @param B The number of boostrap samples used to compute the confidence bands.
#' @param x.eval.grid Input evaluation grid. Defaults to the observed inputs.
#' @param bootstrap.method Bootstrap method to be implemented. See description 
#' for suggested use. Defaults to "nonparametric".
#' @param alpha Specifies the width of the 1-alpha pointwise confidence bands.
#' @param full.ensemble Return the full bootstrap ensemble as an (n x B) matrix.
#' Defaults to \code{FALSE}.
#' @param max_iter Maximum iterations allowed for the trend filtering 
#' ADMM optimization 
#' [\href{http://www.stat.cmu.edu/~ryantibs/papers/fasttf.pdf}{Ramdas & Tibshirani (2015)}]. 
#' Consider increasing this if the trend filtering estimate does not appear to 
#' have fully converged to a reasonable estimate of the signal.
#' @param obj_tol The tolerance used in the ADMM optimization stopping criterion; 
#' when the relative change in the objective function is less than this value, 
#' the algorithm terminates.
#' @param mc.cores Multi-core computing (for speedups): The number of cores to use.
#' Defaults to the number detected on the machine minus 2.
#' @return A list with the following elements:
#' \item{bootstrap.lower.perc.intervals}{Vector of lower bounds for the 1-alpha pointwise confidence band.}
#' \item{bootstrap.upper.perc.intervals}{Vector of upper bounds for the 1-alpha pointwise confidence band.}
#' \item{tf.boot.ensemble}{(Optional) The full bootstrap ensemble as an (n x B) matrix.}
#' @export bootstrap.trendfilter
#' @author Collin A. Politsch, \email{collinpolitsch@@gmail.com}
#' @seealso \code{\link{bootstrap.trendfilter}}
#' @examples 
#' install.packages("devtools")
#' devtools::install_github("statsmaths/glmgen", subdir="R_pkg/glmgen")
#' 
#' # Quasar spectrum example
#' data(quasar)
#' 
#' # SDSS spectra are equally spaced in log10 wavelength space with a separation of 10e-4
#' # Reading in a spectrum file and retrieving the piece of the spectrum in the Lyman-alpha
#' # forest region
#' log.wavelength.scaled <- quasar.spec$col[[2]] * 1000
#' flux <- quasar.spec$col[[1]]
#' wts <- quasar.spec$col[[3]]
#' lya.rest.wavelength <- 1215.67
#' inds <- which(( 10 ^ (log.wavelength.scaled / 1000) ) / (2.953 + 1) < lya.rest.wavelength + 40)
#' log.wavelength.scaled <- log.wavelength.scaled[inds]
#' flux <- flux[inds]
#' wts <- wts[inds]
#' 
#' # Compute SURE loss curve and optimal lambda
#' lambda.grid <- exp(seq(-10,7,length=250))
#' SURE.out <- SURE.trendfilter(log.wavelength.scaled, flux, wts, lambda.grid)
#' lambda.opt <- SURE.out$lambda[which.min(SURE.out$SURE.loss)]
#' 
#' # Fit optimized model
#' fit <- glmgen::trendfilter(log.wavelength.scaled, flux, wts, k = 2, lambda = lambda.opt)
#' 
#' # Plot results
#' wavelength <- 10 ^ (log.wavelength.scaled / 1000)
#' plot(wavelength, flux, type = "l")
#' lines(wavelength, fit$beta, col = "orange", lwd = 2.5)
#' 
#' boot.out <- bootstrap.trendfilter(log.wavelength.scaled, flux, lambda.opt, sigma = sqrt(1/wts),
#'                                   bootstrap.method = "parametric")
#' lines(wavelength, boot.out$bootstrap.lower.perc.intervals, col = "orange", lty = 2, lwd = 2)
#' lines(wavelength, boot.out$bootstrap.upper.perc.intervals, col = "orange", lty = 2, lwd = 2)
#' legend(x = "topleft", lty = c(1,2), col = "orange", lwd = 2, 
#'        legend = c("Trend filtering estimate", "95 percent variability band"))

bootstrap.trendfilter <- function(x, 
                                  y, 
                                  sigma = NULL,
                                  lambda.min,  
                                  k = 2, 
                                  x.eval.grid = x, 
                                  bootstrap.method = "nonparametric", 
                                  alpha = 0.05, 
                                  B = 1000, 
                                  full.ensemble = FALSE,
                                  max_iter = 250, 
                                  obj_tol = 1e-07, 
                                  mc.cores = parallel::detectCores() - 2){
  
  if ( is.null(y) ) stop("y must be specified.")
  if ( is.null(lambda.min) ) stop("lambda.min must be specified.")
  if ( !is.null(x) & (length(x) != length(y)) ) stop("x and y must have same length.")
  if ( !(length(sigma) %in% c(1,length(y),NULL)) ) stop("sigma must either be NULL, a scalar, or the same length as y.")
  
  if ( is.null(sigma) ){
    data <- data.frame(x=x,y=y,wts=1)
  }else{
    data <- data.frame(x=x,y=y,wts=1/sigma^2)
  }
  
  if ( bootstrap.method == "nonparametric" ){
    if ( mc.cores == 1 ){
      tf.boot.ensemble <- matrix(unlist(replicate(B,
                                                  tf.estimator(nonparametric.resampler(data), 
                                                               lambda.opt, 
                                                               k,
                                                               x.eval.grid, 
                                                               max_iter = max_iter, 
                                                               obj_tol = obj_tol
                                                               )
                                                  )
                                        ),
                               ncol = B
                               )
    }else{
      par.func <- function(b){
        boot.tf.estimate <- tf.estimator(nonparametric.resampler(data), 
                                         lambda.opt, 
                                         k, 
                                         x.eval.grid, 
                                         max_iter = max_iter, 
                                         obj_tol = obj_tol
                                         )
        return(boot.tf.estimate)
      }
      tf.boot.ensemble <- matrix(unlist(parallel::mclapply(1:B, par.func, mc.cores = mc.cores)), ncol = B)
    }
  }
  
  if ( bootstrap.method == "parametric" ){
    tf.estimate <- tf.estimator(data, lambda.opt, k, x.eval.grid)
    data$tf.estimate <- tf.estimate
    if ( mc.cores == 1 ){
      tf.boot.ensemble <- matrix(unlist(replicate(B,
                                                  tf.estimator(parametric.sampler(data), 
                                                               lambda.opt, 
                                                               k,
                                                               x.eval.grid, 
                                                               max_iter = max_iter, 
                                                               obj_tol = obj_tol
                                                               )
                                                  )
                                        ),
                                 ncol = B
                                 )
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
      tf.boot.ensemble <- matrix(unlist(replicate(B,
                                                  tf.estimator(wild.sampler(data), 
                                                               lambda.opt, 
                                                               k,
                                                               x.eval.grid, 
                                                               max_iter = max_iter, 
                                                               obj_tol = obj_tol
                                                               )
                                                  )
                                        ),
                               ncol = B
                               )
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
  
  if ( !full.ensemble ){
    return(list(bootstrap.lower.perc.intervals=bootstrap.lower.perc.intervals,
                bootstrap.upper.perc.intervals=bootstrap.upper.perc.intervals))
  }else{
    return(list(bootstrap.lower.perc.intervals=bootstrap.lower.perc.intervals,
                bootstrap.upper.perc.intervals=bootstrap.upper.perc.intervals,
                tf.boot.ensemble = tf.boot.ensemble))
  }
}
