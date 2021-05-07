#' Optimize the trend filtering hyperparameter (with respect to Stein's unbiased 
#' risk estimate)
#'
#' @description \code{SURE.trendfilter} estimates the fixed-input squared error 
#' of a trend filtering estimator on a grid of hyperparameter values via Stein's 
#' unbiased risk estimate [\href{https://projecteuclid.org/journals/annals-of-statistics/volume-9/issue-6/Estimation-of-the-Mean-of-a-Multivariate-Normal-Distribution/10.1214/aos/1176345632.full}{Stein (1981)};
#' \href{http://www.stat.cmu.edu/~larry/=sml/stein.pdf}{Tibshirani & Wasserman 
#' (2015);} 
#' \href{https://academic.oup.com/mnras/article/492/3/4005/5704413}{Politsch et al. (2020)}].
#' @param x The vector of the observed inputs.
#' @param y The vector of the observed outputs.
#' @param weights A vector of weights for the observed outputs. These are
#' defined as \code{weights = 1 / sigma^2}, where \code{sigma} is a vector of 
#' standard errors of the uncertainty in the measured outputs. \code{weights}
#' should either have length equal to 1 (i.e. equiweighted/homoskedastic outputs) 
#' or length equal to \code{length(y)} (i.e. heteroskedastic outputs).
#' @param k (Integer) The degree of the trend filtering estimator. Defaults to 
#' \code{k=2} (quadratic trend filtering).
#' @param lambda A vector of trend filtering hyperparameter values to run the 
#' grid search over. Usually should let them be equally-spaced in log-space (see 
#' Examples). 
#' @param x.eval Grid of inputs to evaluate the optimized trend filtering
#' estimate on. If \code{NULL}, a fine equally-spaced grid is constructed.
#' @param thinning (logical) If \code{TRUE}, then the data are preprocessed so 
#' that a smaller, better conditioned data set is used for fitting. When set to
#' \code{NULL}, the default, function will auto detect whether thinning should 
#' be applied (i.e., cases in which the numerical fitting algorithm will 
#' struggle to converge).
#' @param max_iter Maximum iterations allowed for the trend filtering 
#' convex optimization 
#' [\href{http://www.stat.cmu.edu/~ryantibs/papers/fasttf.pdf}{Ramdas & Tibshirani (2015)}]. 
#' Defaults to \code{max_iter = 200}. Consider increasing this if the trend 
#' filtering estimate does not appear to have fully converged to a reasonable 
#' estimate of the signal.
#' @param obj_tol The tolerance used in the convex optimization stopping 
#' criterion; when the relative change in the objective function is less than 
#' this value, the algorithm terminates. Consider decreasing this if the trend 
#' filtering estimate does not appear to have fully converged to a reasonable 
#' estimate of the signal.
#' @return An object of class 'SURE.trendfilter'. This is a list with the 
#' following elements:
#' \item{x.eval}{The grid of inputs the optimized trend filtering estimate was 
#' evaluated on.}
#' \item{tf.estimate}{The trend filtering estimate of the signal, evaluated on 
#' \code{x.eval}.}
#' \item{validation.method}{"SURE"}
#' \item{lambda}{Vector of hyperparameter values tested during validation.}
#' \item{error}{Vector of estimated SURE errors for the hyperparameter values.}
#' \item{lambda.min}{Hyperparameter value that minimizes the SURE error curve.}
#' \item{df}{Vector of effective degrees of freedom for trend filtering
#' estimators fit during validation.}
#' \item{df.min}{The effective degrees of freedom of the optimally-tuned trend 
#' filtering estimator.}
#' \item{i.min}{The index of \code{lambda} that minimizes the SURE error.}
#' \item{x}{The vector of the observed inputs.}
#' \item{y}{The vector of the observed outputs.}
#' \item{weights}{A vector of weights for the observed outputs. These are
#' defined as \code{weights = 1 / sigma^2}, where \code{sigma} is a vector of 
#' standard errors of the uncertainty in the measured outputs.}
#' \item{fitted.values}{The trend filtering estimate of the signal, evaluated at
#' the observed inputs \code{x}.}
#' \item{residuals}{\code{residuals = y - fitted.values}.}
#' \item{k}{(Integer) The degree of the trend filtering estimator.}
#' \item{thinning}{(logical) If \code{TRUE}, then the data are 
#' preprocessed so that a smaller, better conditioned data set is used for 
#' fitting.}
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
#' @seealso \code{\link{cv.trendfilter}}, \code{\link{bootstrap.trendfilter}}
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
#' ##                    Quasar Lyman-alpha forest example                    ##
#' #############################################################################
#' ##  SDSS spectra are equally spaced in log base-10 wavelength space with a ##
#' ##  separation of 1e-4 log-Angstroms. Given the default trend filtering    ##
#' ##  optimization parameters, it is safer to scale up the inputs in such a  ##
#' ##  scenario. For example, here we scale to unit spacing.                  ##
#' #############################################################################
#' 
#' # Load Lyman-alpha forest spectral observations of an SDSS quasar at redshift 
#' # z = 2.953
#' 
#' data(quasar_spec)
#' data(plotting_utilities)
#' 
#' 
#' # Run the SURE optimization for a quadratic trend filtering estimator, i.e. 
#' # k = 2 (default)
#' 
#' lambda.grid <- exp(seq(-10, 5, length = 150))
#' SURE.obj <- SURE.trendfilter(x = log10.wavelength.scaled, 
#'                              y = flux, 
#'                              weights = weights, 
#'                              lambda = lambda.grid)
#' 
#'                                           
#' # Extract the optimized trend filtering estimate on a fine equally-spaced
#' # grid from the 'SURE.trendfilter' output
#' 
#' lambda.min <- SURE.obj$lambda.min
#' SURE.error <- SURE.obj$error
#' x.eval <- SURE.obj$x.eval                      
#' tf.estimate <- SURE.obj$tf.estimate
#' 
#' 
#' # Transform back to wavelength space
#' 
#' wavelength <- 10 ^ (log10.wavelength.scaled / scale.factor)
#' wavelength.eval <- 10 ^ (x.eval / scale.factor)
#' 
#' 
#' # Plot the results
#'
#' par(mfrow = c(2,1), mar = c(5,4,2.5,1) + 0.1)
#' plot(log(lambda.grid), SURE.error,
#'      main = "SURE error curve", 
#'      xlab = "log(lambda)", ylab = "SURE error")
#' abline(v = log(lambda.min), col = "blue3", lty = 2)
#' text(x = log(lambda.min), y = par("usr")[4], 
#'      labels = "optimal hyperparameter", pos = 1, col = "blue3")
#' 
#' plot(wavelength, flux, type = "l", 
#'      main = "Quasar Lyman-alpha forest", 
#'      xlab = "Observed wavelength (angstroms)", ylab = "flux")
#' lines(wavelength.eval, tf.estimate, col = "orange", lwd = 2.5)
#' legend(x = "topleft", lwd = c(1,2), lty = 1, 
#'        col = c("black","orange"), 
#'        legend = c("Noisy quasar spectrum",
#'                   "Trend filtering estimate"))

SURE.trendfilter <- function(x, 
                             y, 
                             weights, 
                             k = 2L, 
                             lambda, 
                             x.eval = NULL,
                             thinning = NULL,
                             max_iter = 200L, 
                             obj_tol = 1e-06
                             )
  {
  
  if ( !is.numeric(x) ) stop("x must be specified.")
  if ( !is.numeric(y) ) stop("y must be specified.")
  if ( !is.numeric(weights) ) stop("weights are needed in order to compute SURE. If estimates are not available, use cv.trendfilter.")
  if ( !(length(weights) %in% c(1,length(y))) ) stop("weights must either be scalar or same length as y.")
  if ( !is.numeric(lambda) ) stop("a vector of hyperparameter values must be specified.")
  
  lambda <- sort(lambda)

  weights <- case_when(
    length(weights) == 1 ~ rep_len(weights, length(y)),
    length(weights) == length(y) ~ weights
  )
  
  out <- trendfilter(x = x, 
                     y = y,
                     weights = weights, 
                     lambda = lambda,
                     k = k, 
                     thinning = thinning,
                     control = trendfilter.control.list(max_iter = max_iter,
                                                        obj_tol = obj_tol
                                                        )
                     )
                             
  if ( length(lambda) == 1 ){
    SURE.error <- mean( (out$beta - y) ^ 2 ) + (2 * mean(1 / weights) / length(x)) * out$df
  }else{
    SURE.error <- colMeans( (out$beta - y) ^ 2 ) + (2 * mean(1 / weights) / length(x)) * out$df
  }
  
  error <- as.numeric(SURE.error)
  lambda.min <- lambda[which.min(error)]
  
  if ( is.null(x.eval) ) x.eval <- seq(min(x), max(x), length = 1500)
  
  tf.estimate <- glmgen:::predict.trendfilter(out, 
                                              lambda = lambda.min, 
                                              x.new = x.eval
                                              ) %>% as.numeric
  
  fitted.values <- glmgen:::predict.trendfilter(out, 
                                                lambda = lambda.min, 
                                                x.new = x
                                                ) %>% as.numeric
  
  obj <- structure(list(x.eval = x.eval,
                        tf.estimate = tf.estimate,
                        validation.method = "SURE",
                        lambda = lambda, 
                        error = error,
                        lambda.min = lambda.min,
                        df = out$df,
                        df.min = out$df[which.min(error)],
                        i.min = as.integer(which.min(error)),
                        x = x,
                        y = y,
                        weights = weights,
                        fitted.values = fitted.values,
                        residuals = y - fitted.values,
                        k = as.integer(k),
                        thinning = thinning,
                        max_iter = as.integer(max_iter),
                        obj_tol = obj_tol
                        ),
                   class = "SURE.trendfilter"
                   )
  
  return(obj)
}
