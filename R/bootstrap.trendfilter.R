#' Bootstrap the optimized trend filtering estimator to obtain pointwise 
#' variability bands
#'
#' @description \code{bootstrap.trendfilter} implements any of three possible 
#' bootstrap algorithms to obtain pointwise variability bands to accompany the 
#' optimized trend filtering estimate of a signal. 
#' @param x A vector of the observed inputs. 
#' @param y A vector of the observed outputs.
#' @param sigma A vector of measurement standard errors for the observed outputs.
#' @param lambda.min The optimally-tuned trend filtering hyperparameter, e.g. by 
#' minimizing SURE (see \code{\link{SURE.trendfilter}}) or cross validation.
#' @param k The degree of the trend filtering estimator. Defaults to \code{k=2}
#' (quadratic trend filtering).
#' @param B The number of bootstrap samples used to compute the variability 
#' bands.
#' @param x.eval.grid Grid of inputs to evaluate the variability bands on. 
#' Defaults to the observed inputs.
#' @param bootstrap.method Bootstrap method to be used. See Details section
#' below for suggested use. Defaults to "nonparametric".
#' @param alpha Specifies the width of the \code{1-alpha} pointwise variability 
#' bands.
#' @param full.ensemble Return the full bootstrap ensemble as an \code{n x B} 
#' matrix. Defaults to \code{FALSE}.
#' @param max_iter Maximum iterations allowed for the trend filtering 
#' convex optimization 
#' \href{https://stat.cmu.edu/~ryantibs/papers/fasttf.pdf}{
#' Ramdas & Tibshirani (2015)}. 
#' Consider increasing this if the bootstrap estimates do not appear to 
#' have fully converged to a reasonable estimate of the signal.
#' @param obj_tol The tolerance used in the convex optimization stopping 
#' criterion; when the relative change in the objective function is less than 
#' this value, the algorithm terminates. Consider decreasing this if the 
#' bootstrap estimates do not appear to have fully converged to a reasonable 
#' estimate of the signal.
#' @param mc.cores Multi-core computing (for speedups): The number of cores to
#' utilize. Defaults to the number detected on the machine minus 2.
#' @return A list with the following elements:
#' \item{bootstrap.lower.perc.intervals}{Vector of lower bounds for the 1-alpha 
#' pointwise variability band.}
#' \item{bootstrap.upper.perc.intervals}{Vector of upper bounds for the 1-alpha 
#' pointwise variability band.}
#' \item{tf.boot.ensemble}{(Optional) The full bootstrap ensemble as an 
#' \code{n x B} matrix.}
#' @export bootstrap.trendfilter
#' @details The bootstrap method should generally be chosen 
#' according to the following criteria: \cr \cr
#' S1. The inputs are irregularly 
#' sampled –> \code{bootstrap.method = "nonparametric"} \cr \cr
#' S2. The inputs are regularly 
#' sampled and the noise distribution is known –> \code{bootstrap.method = "parametric"} \cr \cr
#' S3. The inputs are regularly sampled and the noise distribution is unknown –> 
#' \code{bootstrap.method = "wild"} \cr
#' @author Collin A. Politsch, \email{collinpolitsch@@gmail.com}
#' @seealso \code{\link{SURE.trendfilter}}
#' @references \enumerate{
#' \item \href{https://academic.oup.com/mnras/article/492/3/4005/5704413}{
#' Politsch et al. (2020). Trend filtering – I. A modern statistical tool for 
#' time-domain astronomy and astronomical spectroscopy} \cr
#' 
#' \item \href{https://academic.oup.com/mnras/article/492/3/4019/5704414}{
#' Politsch et al. (2020). Trend filtering – II. Denoising astronomical signals 
#' with varying degrees of smoothness} \cr
#' 
#' \item \href{https://projecteuclid.org/journals/annals-of-statistics/volume-7/issue-1/Bootstrap-Methods-Another-Look-at-the-Jackknife/10.1214/aos/1176344552.full}{
#' Efron (1979). Bootstrap Methods: Another Look at the Jackknife} \cr
#' 
#' \item \href{https://academic.oup.com/mnras/article/492/3/4019/5704414}{
#' Efron and Tibshirani (1986). Bootstrap Methods for Standard Errors, 
#' Confidence Intervals, and Other Measures of Statistical Accuracy} \cr
#' 
#' \item \href{https://projecteuclid.org/journals/annals-of-statistics/volume-14/issue-4/Jackknife-Bootstrap-and-Other-Resampling-Methods-in-Regression-Analysis/10.1214/aos/1176350142.full}{
#' Wu (1986). Jackknife, Bootstrap and Other Resampling Methods in Regression 
#' Analysis} \cr
#' }
#' @examples 
#' # Quasar spectrum example
#' ## SDSS spectra are equally spaced in log base 10 wavelength space with a 
#' ## separation of 10e-4 logarithmic Angstroms. 
#' 
#' data(quasar_spec)
#' 
#' 
#' # Read in a spectrum of a quasar at redshift z = 2.953 and extract the Lyman-alpha forest.
#' 
#' log.wavelength.scaled <- quasar_spec$col[[2]] * 1000
#' flux <- quasar_spec$col[[1]]
#' wts <- quasar_spec$col[[3]]
#' lya.rest.wavelength <- 1215.67
#' quasar.redshift <- 2.953
#' inds <- which((10^(log.wavelength.scaled/1000))/(quasar.redshift + 1) < lya.rest.wavelength + 40)
#' log.wavelength.scaled <- log.wavelength.scaled[inds]
#' flux <- flux[inds]
#' wts <- wts[inds]
#'
#'
#' # Compute the SURE error curve and the optimal hyperparameter value
#' 
#' lambda.grid <- exp(seq(-10, 7, length = 250))
#' SURE.out <- SURE.trendfilter(x = log.wavelength.scaled, 
#'                              y = flux, 
#'                              sigma = 1 / sqrt(wts), 
#'                              lambda = lambda.grid
#'                              )
#' lambda.min <- SURE.out$lambda.min
#' 
#' 
#' # Fit the optimized trend filtering estimate
#' 
#' fit <- glmgen::trendfilter(x = log.wavelength.scaled, 
#'                            y = flux, 
#'                            weights = wts, 
#'                            k = 2, 
#'                            lambda = lambda.min
#'                            )
#' 
#' 
#' # Plot the results
#' 
#' par(mfrow=c(2,1))
#' plot(log(lambda.grid), SURE.out$SURE.error, xlab = "log(lambda)", ylab = "SURE")
#' abline(v = log(lambda.min), col = "red")
#' wavelength <- 10 ^ (log.wavelength.scaled / 1000)
#' plot(wavelength, flux, type = "l")
#' lines(wavelength, fit$beta, col = "orange", lwd = 2.5)
#' 
#' boot.out <- bootstrap.trendfilter(x = log.wavelength.scaled, 
#'                                   y = flux, 
#'                                   sigma = sqrt(1/wts),
#'                                   lambda.min = lambda.min, 
#'                                   bootstrap.method = "parametric"
#'                                   )
#'                                   
#' lines(wavelength, boot.out$bootstrap.lower.perc.intervals, col = "orange", lty = 2, lwd = 2)
#' lines(wavelength, boot.out$bootstrap.upper.perc.intervals, col = "orange", lty = 2, lwd = 2)
#' legend(x = "topleft", lty = c(1,2), col = "orange", lwd = 2, 
#'        legend = c("Trend filtering estimate", "95 percent variability band"))

#' @importFrom stats quantile
bootstrap.trendfilter <- function(x, 
                                  y, 
                                  sigma,
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
  
  if ( is.null(x) ) stop("x must be specified.")
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
                                                               lambda.min, 
                                                               k,
                                                               edf = NULL,
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
                                         lambda.min, 
                                         k, 
                                         edf = NULL,
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
    tf.estimate <- tf.estimator(data, lambda.min, k, edf = NULL, x.eval.grid)
    data$tf.estimate <- tf.estimate
    if ( mc.cores == 1 ){
      tf.boot.ensemble <- matrix(unlist(replicate(B,
                                                  tf.estimator(parametric.sampler(data), 
                                                               lambda.min, 
                                                               k,
                                                               edf = NULL,
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
        boot.tf.estimate <- tf.estimator(parametric.sampler(data), lambda.min, k, edf = NULL, 
                                         x.eval.grid, max_iter = max_iter, obj_tol = obj_tol)
        return(boot.tf.estimate)
      }
      tf.boot.ensemble <- matrix(unlist(parallel::mclapply(1:B, par.func, mc.cores = mc.cores)), ncol = B)
    }
  }
  
  if ( bootstrap.method == "wild" ){
    data$tf.estimate <- tf.estimator(data, lambda.min, k, edf = NULL, x.eval.grid)
    data$tf.residuals <- data$y - data$tf.estimate
    
    if ( mc.cores == 1 ){
      tf.boot.ensemble <- matrix(unlist(replicate(B,
                                                  tf.estimator(wild.sampler(data), 
                                                               lambda.min, 
                                                               k,
                                                               edf = NULL,
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
        boot.tf.estimate <- tf.estimator(wild.sampler(data), lambda.min, k, edf = NULL,
                                         x.eval.grid, max_iter = max_iter, obj_tol = obj_tol)
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


tf.estimator <- function(data, lambda, k, edf = NULL, x.eval.grid, max_iter = 250, obj_tol = 1e-06){
  tf.fit <- glmgen::trendfilter(data$x, 
                                data$y, 
                                data$wts, 
                                k = 2, 
                                lambda = lambda, 
                                control = glmgen::trendfilter.control.list(max_iter = max_iter, 
                                                                           obj_tol = obj_tol
                                )
  )
  tf.estimate <- as.numeric(suppressWarnings(glmgen:::predict.trendfilter(tf.fit, x.new = x.eval.grid, lambda = lambda)))
  return(tf.estimate)
}


#' @importFrom dplyr slice_sample
nonparametric.resampler <- function(data){
  resampled.data <- slice_sample(data, size = nrow(data), replace = TRUE)
  return(resampled.data)
}


#' @importFrom stats rnorm
parametric.sampler <- function(data){
  boot.sample <- data$tf.estimate + rnorm(nrow(data), sd = 1 / sqrt(data$wts))
  return(data.frame(x=data$x,y=boot.sample,wts=data$wts))
}


wild.sampler <- function(data){
  wild.boot.resids <- data$tf.residuals * sample(x = c((1+sqrt(5))/2, 1-sqrt(5)/2), size = nrow(data), replace = T, 
                                                 prob = c((1+sqrt(5))/(2*sqrt(5)),(sqrt(5)-1)/(2*sqrt(5))))
  wild.boot.sample <- data$tf.estimate + wild.boot.resids
  return(data.frame(x=data$x,y=wild.boot.sample,wts=data$wts))
}
