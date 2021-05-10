#' Optimize the trend filtering hyperparameter (with respect to Stein's unbiased 
#' risk estimate)
#'
#' @description \code{SURE.trendfilter} estimates the fixed-input squared error 
#' of a trend filtering estimator (via Stein's unbiased risk estimate) over a 
#' grid of hyperparameter values and returns the optimized estimator.
#' @param x The vector of observed values of the input variable (a.k.a. the 
#' predictor, covariate, explanatory variable, regressor, independent variable, 
#' control variable, etc.)
#' @param y The vector of observed values of the output variable (a.k.a. the
#' response, target, outcome, regressand, dependent variable, etc.).
#' @param weights (Must be passed in order to optimize with respect to SURE):
#' A vector of weights for the observed outputs. These are defined as 
#' \code{weights = 1 / sigma^2}, where \code{sigma} is a vector of standard 
#' errors of the uncertainty in the output measurements. \code{weights} should 
#' either have length equal to 1 (corresponding to observations with a constant 
#' (scalar) variance of \code{sigma = 1/sqrt(weights)}) or length equal to 
#' \code{length(y)} (i.e. heteroskedastic outputs). If reliable estimates of
#' the error variance are not available, use \code{cv.trendfilter} instead.
#' @param k The degree of the trend filtering estimator. Defaults to 
#' \code{k=2} (quadratic trend filtering). Must be one of \code{k = 0,1,2,3},
#' although \code{k=3} is discouraged due to algorithmic instability (and is
#' visually indistinguishable from \code{k=2} anyway).
#' @param nlambda The number of trend filtering hyperparameter values 
#' to run the grid search over. If \code{lambda = NULL}, a grid of length 
#' \code{nlambda} is constructed internally.
#' @param lambda Overrides \code{nlambda} if passed. A user-supplied vector of 
#' trend filtering hyperparameter values to run the grid search over. Usually, 
#' let them be equally-spaced in log-space (see Examples), and good to 
#' provide them in descending order.
#' @param n.eval The length of the equally-spaced \code{x} grid to evaluate the 
#' optimized trend filtering estimate on.
#' @param x.eval Overrides \code{n.eval} if passed. A user-supplied grid of 
#' inputs to evaluate the optimized trend filtering estimate on. 
#' @param thinning If \code{TRUE}, then the data are preprocessed so 
#' that a smaller, better conditioned data set is used for fitting. When left
#' \code{NULL}, the default, the optimization will automatically detect whether 
#' thinning should be applied (i.e., cases in which the numerical fitting 
#' algorithm will struggle to converge).
#' @param max_iter Maximum iterations allowed for the trend filtering 
#' convex optimization 
#' [\href{http://www.stat.cmu.edu/~ryantibs/papers/fasttf.pdf}{Ramdas & Tibshirani (2015)}]. 
#' Defaults to \code{max_iter = 250L}. Consider increasing this if the trend 
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
#' \item{tf.estimate}{The optimized trend filtering estimate of the signal, 
#' evaluated on \code{x.eval}.}
#' \item{validation.method}{"SURE"}
#' \item{lambda}{Vector of hyperparameter values tested during validation
#' (always returned in descending order).}
#' \item{error}{Vector of SURE error estimates corresponding to the *descending* 
#' set of lambda values tested during validation.}
#' \item{lambda.min}{Hyperparameter value that minimizes the SURE error curve.}
#' \item{df}{Vector of effective degrees of freedom for all trend filtering
#' estimators fit during validation.}
#' \item{df.min}{The effective degrees of freedom of the optimally-tuned trend 
#' filtering estimator.}
#' \item{i.min}{The index of \code{lambda} (descending order) that minimizes 
#' the SURE error curve.}
#' \item{x}{The vector of the observed inputs.}
#' \item{y}{The vector of the observed outputs.}
#' \item{weights}{A vector of weights for the observed outputs. These are
#' defined as \code{weights = 1 / sigma^2}, where \code{sigma} is a vector of 
#' standard errors of the uncertainty in the measured outputs.}
#' \item{fitted.values}{The optimized trend filtering estimate of the signal, 
#' evaluated at the observed inputs \code{x}.}
#' \item{residuals}{\code{residuals = y - fitted.values}.}
#' \item{k}{The degree of the trend filtering estimator.}
#' \item{thinning}{If \code{TRUE}, then the data are preprocessed so that a 
#' smaller, better conditioned data set is used for fitting.}
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
#' 
#' # Load Lyman-alpha forest spectral observations of an SDSS quasar at redshift 
#' # z = 2.953. SDSS spectra are equally spaced in log10 wavelength space.
#' 
#' data(quasar_spec)
#' data(plotting_utilities)
#' 
#' 
#' # Run the SURE optimization for a quadratic trend filtering estimator, i.e. 
#' # k = 2 (default)
#' 
#' SURE.obj <- SURE.trendfilter(x = log10.wavelength, 
#'                              y = flux, 
#'                              weights = weights)
#' 
#' 
#' # Extract the SURE error curve and optimized trend filtering estimate from 
#' # the output
#' 
#' log.lambda <- log(SURE.obj$lambda)
#' error <- SURE.obj$error
#' log.lambda.min <- log(SURE.obj$lambda.min)
#' 
#' log10.wavelength.eval <- SURE.obj$x.eval
#' tf.estimate <- SURE.obj$tf.estimate
#' 
#' 
#' # Transform the inputs to wavelength space (in Angstroms)
#' 
#' wavelength <- 10 ^ (log10.wavelength)
#' wavelength.eval <- 10 ^ (log10.wavelength.eval)
#' 
#' 
#' # Plot the results
#'
#' par(mfrow = c(2,1), mar = c(5,4,2.5,1) + 0.1)
#' plot(x = log.lambda, y = error, main = "SURE error curve", 
#'      xlab = "log(lambda)", ylab = "SURE error")
#' abline(v = log.lambda.min, lty = 2, col = "blue3")
#' text(x = log.lambda.min, y = par("usr")[4], 
#'      labels = "optimal hyperparameter", pos = 1, col = "blue3")
#' 
#' plot(x = wavelength, y = flux, type = "l", main = "Quasar Lyman-alpha forest", 
#'      xlab = "Observed wavelength (Angstroms)", ylab = "Flux")
#' lines(wavelength.eval, tf.estimate, col = "orange", lwd = 2.5)
#' legend(x = "topleft", lwd = c(1,2), lty = 1, col = c("black","orange"), 
#'        legend = c("Noisy quasar Lyman-alpha forest", "Trend filtering estimate"))

SURE.trendfilter <- function(x, 
                             y, 
                             weights, 
                             k = 2L, 
                             nlambda = 250L, 
                             lambda = NULL,
                             n.eval = 1500L,
                             x.eval = NULL,
                             thinning = NULL,
                             max_iter = 500L, 
                             obj_tol = 1e-10
                             )
  {
  
  if ( missing(x) || is.null(x) ) stop("x must be passed.")
  if ( missing(y) || is.null(y) ) stop("y must be passed.")
  if ( length(x) != length(y) ) stop("x and y must have the same length.")
  if ( missing(weights) | !is.numeric(weights) ){
    stop("weights are needed in order to compute SURE. If estimates are not available, use cv.trendfilter.")
  }
  if ( any(weights == 0L) ) stop("cannot pass zero weights.")
  if ( !(length(weights) %in% c(1,length(y))) ) stop("weights must either be have length 1 or length(y).")
  if ( length(y) < k + 2 ) stop("y must have length >= k+2 for kth order trend filtering.")
  if ( k < 0 || k != round(k) ) stop("k must be a nonnegative integer. k=2 recommended")
  if ( k > 3 ) stop("Large k leads to generally worse conditioning; k=0,1,2 are the most stable choices.")
  if ( k == 3 ) warning("k=3 can have poor conditioning; k=2 is more stable and visually identical.")
  if ( is.null(lambda) ){
    if ( nlambda != round(nlambda) || nlambda < 25L ) stop("nlambda must be a positive integer >= 25.")
  }else{
    if ( min(lambda) < 0L ) stop("All specified lambda values must be nonnegative.")
    if ( length(lambda) < 25L ) stop("lambda must be have length >= 25.")
  }
  if ( length(weights) == 1 ) weights <- rep(weights, length(x))
  
  data <- tibble(x, y, weights) %>% 
    arrange(x) %>% 
    drop_na

  n.obs <- nrow(data)
  x.scale <- mean(diff(data$x))
  y.scale <- mean(data$y) / 10
  x <- data$x / x.scale
  y <- data$y / y.scale
  weights <- y.scale ^ 2 * data$weights
  
  if ( is.null(lambda) ){
    lambda <- seq(10, -10, length = nlambda) %>% exp 
  }else{
    lambda <- sort(lambda, decreasing = TRUE)
  }

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
  
  if ( length(lambda) == 1 ){
    SURE.error <- y.scale ^ 2 * mean( (out$beta - y) ^ 2 ) + 2 * out$df / n.obs * mean(y.scale ^ 2 / weights)
  }else{
    SURE.error <- y.scale ^ 2 * colMeans( (out$beta - y) ^ 2 ) + 2 * out$df / n.obs * mean(y.scale ^ 2 / weights)
  }

  error <- as.numeric(SURE.error)
  i.min <- as.integer(which.min(error))
  lambda.min <- lambda[i.min]
  
  if ( is.null(x.eval) ){
    x.eval <- seq(min(x), max(x), length = n.eval)
  }else{
    x.eval <- sort(x.eval) / x.scale
  }
  
  tf.estimate <- glmgen:::predict.trendfilter(out, 
                                              lambda = lambda.min, 
                                              x.new = x.eval
                                              ) %>% as.numeric
  
  fitted.values <- glmgen:::predict.trendfilter(out, 
                                                lambda = lambda.min, 
                                                x.new = x
                                                ) %>% as.numeric
  
  obj <- structure(list(x.eval = x.eval * x.scale,
                        tf.estimate = tf.estimate * y.scale,
                        validation.method = "SURE",
                        lambda = lambda, 
                        error = error,
                        lambda.min = lambda.min,
                        df = out$df,
                        df.min = out$df[i.min],
                        i.min = i.min,
                        x = x * x.scale,
                        y = y * y.scale,
                        weights = weights / y.scale ^ 2,
                        fitted.values = fitted.values * y.scale,
                        residuals = (y - fitted.values) * y.scale,
                        k = as.integer(k),
                        thinning = thinning,
                        max_iter = as.integer(max_iter),
                        obj_tol = obj_tol
                        ),
                   class = "SURE.trendfilter"
                   )
  
  return(obj)
}
