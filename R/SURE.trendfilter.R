#' Optimize the trend filtering hyperparameter (with respect to Stein's unbiased 
#' risk estimate)
#'
#' @description \code{SURE.trendfilter} computes the Stein's unbiased risk 
#' estimate of fixed-input mean-squared prediction error on a grid of 
#' hyperparameter values and returns the full error curve and the optimal
#' hyperparameter value.
#' @param x A vector of the observed inputs. If \code{NULL}, then we assume
#' unit spacing.
#' @param y A vector of the observed outputs.
#' @param sigma A vector of measurement standard errors for the observed outputs.
#' @param lambda A vector of trend filtering hyperparameter values to run the 
#' grid search over.
#' @param k The degree of the trend filtering estimator. Defaults to \code{k=2}
#' (quadratic trend filtering).
#' @param max_iter Maximum iterations allowed for the trend filtering 
#' ADMM optimization 
#' [\href{http://www.stat.cmu.edu/~ryantibs/papers/fasttf.pdf}{Ramdas & Tibshirani (2015)}]. 
#' Consider increasing this if the trend filtering estimate does not appear to 
#' have fully converged to a reasonable estimate of the signal.
#' @param obj_tol The tolerance used in the ADMM optimization stopping criterion; 
#' when the relative change in the objective function is less than this value, 
#' the algorithm terminates.
#' @return A list with the following elements:
#' \item{lambda}{Vector of hyperparameter values tested.}
#' \item{SURE.loss}{Vector of estimated SURE errors for hyperparameter values.}
#' \item{lambda.min}{Hyperparameter value that minimizes the SURE curve.}
#' @export SURE.trendfilter
#' @author Collin A. Politsch, \email{collinpolitsch@@gmail.com}
#' @seealso \link{bootstrap.trendfilter}
#' @examples 
#' # Quasar spectrum example
#' # SDSS spectra are equally spaced in log base 10 wavelength space with a 
#' # separation of 10e-4 logarithmic Angstroms. 
#' 
#' data(quasar)
#' 
#' 
#' # Read in a spectrum of a quasar at redshift z = 2.953 and extract the
#' # Lyman-alpha forest.
#' 
#' log.wavelength.scaled <- quasar.spec$col[[2]] * 1000
#' flux <- quasar.spec$col[[1]]
#' wts <- quasar.spec$col[[3]]
#' lya.rest.wavelength <- 1215.67
#' inds <- which((10^(log.wavelength.scaled/1000))/(2.953 + 1) < lya.rest.wavelength + 40)
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
#' plot(log(lambda.grid), SURE.out$SURE.loss, xlab = "log(lambda)", ylab = "SURE")
#' abline(v = log(lambda.min), col = "red")
#' wavelength <- 10 ^ (log.wavelength.scaled / 1000)
#' plot(wavelength, flux, type = "l")
#' lines(wavelength, fit$beta, col = "orange", lwd = 2.5)

SURE.trendfilter <- function(x, 
                             y, 
                             sigma, 
                             lambda, 
                             k = 2, 
                             max_iter = 250, 
                             obj_tol = 1e-06
                             ){
  
  if ( is.null(x) ) stop("x must be specified.")
  if ( is.null(y) ) stop("y must be specified.")
  if ( is.null(sigma) ) stop("sigma must be provided in order to compute SURE.")
  if ( is.null(lambda) ) stop("lambda must be specified.")
  if ( !(length(sigma) %in% c(1,length(y))) ) stop("sigma must either be scalar or same length as y.")
  
  sigma <- ifelse(length(sigma) == 1, rep(sigma, length(y)), sigma)
  wts <- 1/sigma^2
  out <- glmgen::trendfilter(x = x, 
                             y = y,
                             weights = wts, 
                             k = k, 
                             lambda = lambda,
                             control = glmgen::trendfilter.control.list(max_iter = max_iter, 
                                                                        obj_tol = obj_tol)
                             )
  if ( length(lambda) == 1 ){
    SURE.loss <- mean( (out$beta - y)^2 ) + (2 * mean(1/wts) / length(x)) * out$df
  }
  if ( length(lambda) > 1 ){
    SURE.loss <- as.numeric(colMeans( (out$beta - y)^2 )) + (2 * mean(1/wts) / length(x)) * out$df
  }
  return(list(lambda = lambda, 
              SURE.loss = SURE.loss, 
              lambda.min = lambda[which.min(SURE.loss)],
              df.min = 
              )
         )
}
