#' Optimize the trend filtering hyperparameter (with respect to Stein's unbiased 
#' risk estimate)
#'
#' @description \code{SURE.trendfilter} computes the Stein's unbiased risk 
#' estimate of fixed-input mean-squared prediction error (MSPE) on a grid of 
#' hyperparameter values.
#' @param x A vector of the observed inputs.
#' @param y A vector of the observed outputs.
#' @param weights A vector of weights for the observed outputs. These are
#' defined as \code{weights = 1 / sigma^2}, where \code{sigma} is a vector of 
#' standard errors of the uncertainty in the measured outputs. \code{weights}
#' should either have length equal to 1 (i.e. equiweighted/homoskedastic outputs) 
#' or length equal to \code{length(y)} (i.e. heteroskedastic outputs).
#' @param k (Integer) The degree of the trend filtering estimator. Defaults to 
#' \code{k = 2} (quadratic trend filtering).
#' @param lambda A vector of trend filtering hyperparameter values to run the 
#' grid search over. Usually, let them be equally-spaced in log-space (see 
#' Examples). 
#' @param max_iter Maximum iterations allowed for the trend filtering 
#' convex optimization 
#' [\href{http://www.stat.cmu.edu/~ryantibs/papers/fasttf.pdf}{Ramdas & Tibshirani (2015)}]. 
#' Defaults to \code{max_iter = 150}. Consider increasing this if the trend 
#' filtering estimate does not appear to have fully converged to a reasonable 
#' estimate of the signal.
#' @param obj_tol The tolerance used in the convex optimization stopping 
#' criterion; when the relative change in the objective function is less than 
#' this value, the algorithm terminates. Consider decreasing this if the trend 
#' filtering estimate does not appear to have fully converged to a reasonable 
#' estimate of the signal.
#' @return An object of class \code{SURE.trendfilter}. This is a list with the 
#' following elements:
#' \item{error}{Vector of estimated SURE errors for hyperparameter values.}
#' \item{lambda}{Vector of hyperparameter values tested.}
#' \item{lambda.min}{Hyperparameter value that minimizes the SURE error curve.}
#' \item{df}{Vector of effective degrees of freedom for all trend filtering
#' estimators with hyperparameters \code{lambda}.}
#' \item{df.min}{The effective degrees of freedom of the optimally-tuned trend 
#' filtering estimator.}
#' \item{i.min}{The index of \code{lambda} that minimizes the SURE error.}
#' \item{x}{A vector of the observed inputs.}
#' \item{y}{A vector of the observed outputs.}
#' \item{weights}{A vector of weights for the observed outputs. These are
#' defined as \code{weights = 1 / sigma^2}, where \code{sigma} is a vector of 
#' standard errors of the uncertainty in the measured outputs. \code{weights}
#' should either have length equal to 1 (i.e. equiweighted/homoskedastic outputs) 
#' or length equal to \code{length(y)} (i.e. heteroskedastic outputs).}
#' \item{k}{(Integer) The degree of the trend filtering estimator.}
#' \item{max_iter}{Maximum iterations allowed for the trend filtering 
#' convex optimization 
#' [\href{http://www.stat.cmu.edu/~ryantibs/papers/fasttf.pdf}{Ramdas & Tibshirani (2015)}]. 
#' Consider increasing this if the trend filtering estimate does not appear to 
#' have fully converged to a reasonable estimate of the signal.}
#' \item{obj_tol}{The tolerance used in the convex optimization stopping 
#' criterion; when the relative change in the objective function is less than 
#' this value, the algorithm terminates. Consider decreasing this if the trend 
#' filtering estimate does not appear to have fully converged to a reasonable 
#' estimate of the signal.}
#' @export SURE.trendfilter
#' @author Collin A. Politsch, \email{collinpolitsch@@gmail.com}
#' @seealso \link{bootstrap.trendfilter}
#' @references \enumerate{
#' \item \href{https://academic.oup.com/mnras/article/492/3/4005/5704413}{
#' Politsch et al. (2020). Trend filtering – I. A modern statistical tool for 
#' time-domain astronomy and astronomical spectroscopy} \cr
#' 
#' \item \href{https://academic.oup.com/mnras/article/492/3/4019/5704414}{
#' Politsch et al. (2020). Trend filtering – II. Denoising astronomical signals 
#' with varying degrees of smoothness} \cr
#' 
#' \item \href{http://www.stat.cmu.edu/~larry/=sml/stein.pdf}{Tibshirani 
#' and Wasserman (Course notes). Stein’s Unbiased Risk Estimate} \cr
#' 
#' \item \href{https://www.tandfonline.com/doi/abs/10.1198/016214504000000692}{Efron 
#' (2004). The Estimation of Prediction Error: Covariance Penalties and Cross-Validation} \cr
#' 
#' \item \href{https://projecteuclid.org/journals/annals-of-statistics/volume-9/issue-6/Estimation-of-the-Mean-of-a-Multivariate-Normal-Distribution/10.1214/aos/1176345632.full}{Stein (1981).
#' Estimation of the Mean of a Multivariate Normal Distribution}
#' }
#' @examples 
#' #############################################################################
#' ##################### Quasar Lyman-alpha forest example #####################
#' #############################################################################
#' 
#' # SDSS spectra are equally spaced in log base-10 wavelength space with a 
#' # separation of 10e-4 log-Angstroms. Given the default trend filtering 
#' # optimization parameters, it is safer to scale up the inputs in such a 
#' # scenario. Here, we scale to unit spacing.
#' 
#' # Read in an SDSS spectrum of a quasar at redshift z = 2.953 and extract the 
#' # Lyman-alpha forest.
#' 
#' data(quasar_spec)
#' lya.rest <- 1215.67
#' quasar.redshift <- 2.953
#' 
#' log.wavelength.scaled <- quasar_spec$col[[2]] * 1000
#' flux <- quasar_spec$col[[1]]
#' weights <- quasar_spec$col[[3]]
#' 
#' inds <- which((10^(quasar_spec$col[[2]]))/(quasar.redshift + 1) < lya.rest)
#' x <- log.wavelength.scaled[inds]
#' y <- flux[inds]
#' weights <- weights[inds]
#'
#'
#' # Run the SURE optimization for a quadratic trend filtering estimator, i.e. 
#' # k = 2 (recommended)
#' 
#' set.seed(1)
#' lambda.grid <- exp(seq(-10, 7, length = 250))
#' SURE.obj <- SURE.trendfilter(x = x, 
#'                              y = y, 
#'                              weights = weights, 
#'                              k = 2,
#'                              lambda = lambda.grid
#'                              )
#' lambda.min <- SURE.obj$lambda.min
#' 
#' 
#' # Fit the optimized trend filtering model and get the estimates on an fine
#' # equally-spaced input grid
#' 
#' model <- trendfilter(x = x,
#'                      y = y, 
#'                      weights = weights,
#'                      k = 2, 
#'                      lambda = lambda.min
#'                      )
#'                      
#' x.eval.grid <- seq(min(x), max(x), length = 1500)
#' tf.estimate <- predict(model, x.new = x.eval.grid)
#' 
#' 
#' # Plot the results
#'
#' par(mfrow = c(2,1), mar = c(5,4,2.5,1) + 0.1)
#' 
#' plot(log(lambda.grid), SURE.obj$error,
#'      main = "SURE error curve", 
#'      xlab = "log(lambda)", ylab = "SURE error")
#' abline(v = log(lambda.min), col = "blue3", lty = 2)
#' text(x = log(lambda.min), y = par("usr")[4], 
#'      labels = "optimal hyperparameter", pos = 1, col = "blue3")
#'      
#'      
#' # Transform back to wavelength space
#' wavelength <- 10 ^ (x / 1000)
#' wavelength.eval.grid <- 10 ^ (x.eval.grid / 1000)
#' 
#' plot(wavelength, y, type = "l", 
#'      main = "Quasar Lyman-alpha forest", 
#'      xlab = "Observed wavelength (angstroms)", ylab = "flux")
#' lines(wavelength.eval.grid, tf.estimate, col = "orange", lwd = 2.5)

#' @import glmgen
SURE.trendfilter <- function(x, 
                             y, 
                             weights, 
                             k = 2L, 
                             lambda, 
                             max_iter = 150L, 
                             obj_tol = 1e-06
                             )
  {
  
  if ( !is.numeric(x) ) stop("x must be specified.")
  if ( !is.numeric(y) ) stop("y must be specified.")
  if ( !is.numeric(weights) ) stop("weights are needed in order to compute SURE. If estimates are not available, use cross validation.")
  if ( !(length(weights) %in% c(1,length(y))) ) stop("weights must either be scalar or same length as y.")
  if ( !is.numeric(lambda) ) stop("lambda must be specified.")
  if ( length(weights) == 1 ) weights <- rep(weights, times = length(y))
  
  out <- trendfilter(x = x, 
                     y = y,
                     weights = weights, 
                     lambda = lambda,
                     k = k, 
                     control = glmgen::trendfilter.control.list(max_iter = max_iter, 
                                                                obj_tol = obj_tol
                                                                )
                     )
                             
  if ( length(lambda) == 1 ){
    SURE.error <- mean( (out$beta - y) ^ 2 ) + (2 * mean(1 / weights) / length(x)) * out$df
  }
  if ( length(lambda) > 1 ){
    SURE.error <- colMeans( (out$beta - y) ^ 2 ) + (2 * mean(1 / weights) / length(x)) * out$df
  }
  
  error <- as.numeric(SURE.error)
  
  out.arg <- structure(list(error = error,
                            lambda = lambda, 
                            lambda.min = lambda[which.min(error)],
                            df = out$df,
                            df.min = out$df[which.min(error)],
                            i.min = which.min(error),
                            x = x,
                            y = y,
                            weights = weights,
                            k = as.integer(k),
                            max_iter = max_iter,
                            obj_tol = obj_tol
                            ),
                       class = "SURE.trendfilter"
                       )
  
  return(out.arg)
}
