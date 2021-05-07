#' Optimize the trend filtering hyperparameter (by V-fold cross validation)
#'
#' @description \code{cv.trendfilter} performs V-fold cross validation to
#' estimate the random-input squared error of a trend filtering estimator on a 
#' grid of hyperparameter values.
#' @param x The vector of the observed inputs.
#' @param y The vector of the observed outputs.
#' @param weights A vector of weights for the observed outputs. These are
#' defined as \code{weights = 1/sigma^2}, where \code{sigma} is a vector of 
#' standard errors of the uncertainty in the measured outputs. Defaults to 
#' \code{NULL}, in which case the data all have unit weighting. If an
#' argument is supplied, it should either have length equal to 1 
#' (implying homoskedastic outputs) or a vector with length equal to 
#' \code{length(y)} (i.e. heteroskedastic outputs).
#' @param k (Integer) The degree of the trend filtering estimator. Defaults to 
#' \code{k=2} (quadratic trend filtering).
#' @param lambda A vector of trend filtering hyperparameter values to run the 
#' grid search over. Usually, let them be equally-spaced in log-space (see 
#' Examples). 
#' @param V The number of folds the data are split into for the V-fold cross
#' validation. Defaults to \code{V=10} (recommended).
#' \code{V=length(x)} is equivalent to leave-one-out cross validation.
#' @param validation.error.type Type of error to optimize during cross
#' validation. One of c("MSE","MAD"), i.e. either mean-squared error or mean 
#' absolute deviations error. Defaults to \code{"MSE"}.
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
#' @param mc.cores Multi-core computing (for speedups): The number of cores to
#' utilize. If 4 or more cores are detected, then the default is to utilize the
#' minimum of \code{available.cores - 2} and \code{V}. Else, \code{mc.cores = 1}.
#' @return An object of class 'cv.trendfilter'. This is a list with the 
#' following elements:
#' \item{x.eval}{The grid of inputs the optimized trend filtering estimate was 
#' evaluated on.}
#' \item{tf.estimate}{The optimizied trend filtering estimate of the signal, 
#' evaluated on \code{x.eval}.}
#' \item{validation.method}{"cv"}
#' \item{V}{The number of folds the data are split into for the V-fold cross
#' validation.}
#' \item{validation.error.type}{Type of validation loss. One of c("MSE","MAD").}
#' \item{lambda}{Vector of hyperparameter values tested during validation.}
#' \item{error}{Vector of cross validation errors for the given hyperparameter 
#' values.}
#' \item{se.error}{The standard errors of the cross validation errors.
#' These are particularly useful for implementing the ``1-standard-error rule''. 
#' The 1-SE rule favors a smoother trend filtering estimate by, instead of 
#' using the hyperparameter that minimizes the CV error, instead uses the 
#' largest hyperparameter that has a CV error within 1 standard error of the
#' smallest CV error.}
#' \item{lambda.min}{Hyperparameter value that minimizes the SURE error curve.}
#' \item{lambda.1se}{The largest hyperparameter value that is still within one
#' standard error of the minimum hyperparameter's cross validation error.}
#' \item{df}{Vector of effective degrees of freedom for trend filtering
#' estimators fit during validation.}
#' \item{df.min}{The effective degrees of freedom of the optimally-tuned trend 
#' filtering estimator.}
#' \item{df.1se}{The effective degrees of freedom of the 1-stand-error rule
#' trend filtering estimator.}
#' \item{i.min}{The index of \code{lambda} that minimizes the cross validation 
#' error.}
#' \item{i.1se}{The index of \code{lambda} that gives the largest hyperparameter
#' value that has a cross validation error within 1 standard error of the 
#' minimum of the cross validation error curves.}
#' \item{x}{The vector of the observed inputs.}
#' \item{y}{The vector of the observed outputs.}
#' \item{weights}{A vector of weights for the observed outputs. These are
#' defined as \code{weights = 1/sigma^2}, where \code{sigma} is a vector of 
#' standard errors of the uncertainty in the measured outputs.}
#' \item{fitted.values}{The trend filtering estimate of the signal, evaluated at
#' the observed inputs \code{x}.}
#' \item{residuals}{\code{residuals = y - fitted.values}.}
#' \item{k}{(Integer) The degree of the trend filtering estimator.}
#' \item{thinning}{(logical) If \code{TRUE}, then the data are preprocessed so 
#' that a smaller, better conditioned data set is used for fitting. When set to
#' \code{NULL}, the default, function will auto detect whether thinning should 
#' be applied (i.e., cases in which the numerical fitting algorithm will 
#' struggle to converge).}
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
#' @export cv.trendfilter
#' @author Collin A. Politsch, \email{collinpolitsch@@gmail.com}
#' @seealso \code{\link{SURE.trendfilter}}, \code{\link{bootstrap.trendfilter}}
#' @references \enumerate{
#' \item \href{https://academic.oup.com/mnras/article/492/3/4005/5704413}{
#' Politsch et al. (2020). Trend filtering – I. A modern statistical tool for 
#' time-domain astronomy and astronomical spectroscopy} \cr
#' 
#' \item \href{https://academic.oup.com/mnras/article/492/3/4019/5704414}{
#' Politsch et al. (2020). Trend filtering – II. Denoising astronomical signals 
#' with varying degrees of smoothness} \cr
#' 
#' \item 
#' 
#' \item 
#' 
#' \item 
#' }
#' @examples 
#' #############################################################################
#' ##                    Quasar Lyman-alpha forest example                    ##
#' #############################################################################
#' ##  SDSS spectra are equally spaced in log base-10 wavelength space with a ##
#' ##  separation of 10e-4 log-Angstroms. Given the default trend filtering   ##
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
#' # Run the cross validation for a quadratic trend filtering estimator, i.e. 
#' # k = 2 (default)
#' 
#' lambda.grid <- exp(seq(-10, 5, length = 150))
#' cv.obj <- cv.trendfilter(x = log10.wavelength.scaled,
#'                          y = flux,
#'                          weights = weights,
#'                          lambda = lambda.grid)
#' 
#'                                           
#' # Extract the optimized trend filtering estimate on a fine equally-spaced
#' # grid from the 'cv.trendfilter' output
#' 
#' lambda.min <- cv.obj$lambda.min
#' cv.error <- cv.obj$error
#' x.eval <- cv.obj$x.eval                      
#' tf.estimate <- cv.obj$tf.estimate
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
#' plot(log(lambda.grid), cv.error,
#'      main = "CV error curve", 
#'      xlab = "log(lambda)", ylab = "CV error")
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

#' @import magrittr
#' @import glmgen
#' @importFrom parallel mclapply
#' @importFrom matrixStats rowSds
cv.trendfilter <- function(x, 
                           y, 
                           weights = NULL, 
                           k = 2L, 
                           lambda, 
                           V = 10L,
                           validation.error.type = c("MSE","MAD"),
                           x.eval = NULL,
                           thinning = NULL,
                           max_iter = 200L, 
                           obj_tol = 1e-06,
                           mc.cores = max(c(parallel::detectCores() - 2), 1)
                           )
  {
  
  if ( !is.numeric(x) ) stop("x must be specified.")
  if ( !is.numeric(y) ) stop("y must be specified.")
  if ( !(length(weights) %in% c(1,length(y))) ) stop("weights must either be scalar or same length as y.")
  if ( !is.numeric(lambda) ) stop("a vector of hyperparameter values must be specified.")
  if ( length(lambda) < 5 ) stop("your hyperparameter grid is too coarse.")
  
  weights <- case_when(
    is.null(weights) ~ rep_len(1, length(y)),
    length(weights) == 1 ~ rep_len(weights, length(y)),
    length(weights) == length(y) ~ weights
  )
  
  if ( is.null(x.eval) ) x.eval <- seq(min(x), max(x), length = 1500)
  lambda <- sort(lambda)
  validation.error.type <- match.arg(validation.error.type)
  mc.cores <- ifelse(mc.cores > V, V, mc.cores)
  
  obj <- structure(list(x.eval = x.eval,
                        validation.method = "cv",
                        V = as.integer(V),
                        validation.error.type = validation.error.type,
                        lambda = lambda, 
                        x = x,
                        y = y,
                        weights = weights,
                        k = as.integer(k),
                        thinning = thinning,
                        max_iter = as.integer(max_iter),
                        obj_tol = obj_tol
                        ),
                   class = "cv.trendfilter"
                   )
  
  data.folded <- tibble(x, y, weights) %>% 
    group_split( sample( rep_len(1:V, length(x)) ), .keep = FALSE )

  cv.out <- mclapply(1:V, 
                     FUN = trendfilter.validate, 
                     data.folded = data.folded, 
                     obj = obj, 
                     mc.cores = mc.cores
                     ) %>%
    unlist %>%
    matrix(ncol = V)
  
  obj <- c(obj, list(error = rowMeans(cv.out),
                     se.error = matrixStats::rowSds(cv.out) / sqrt(V)
                     )
           )
 
  obj$i.min <- which.min(obj$error)
  obj$lambda.min <- obj$lambda[obj$i.min]
  obj$i.1se <- obj %$% which(error <= error[i.min] + se.error[i.min]) %>% max
  obj$lambda.1se <- obj$lambda[obj$i.1se]

  out <- obj %$%
    trendfilter(x = x, 
                y = y,
                weights = weights,
                lambda = obj$lambda,
                k = obj$k, 
                thinning = obj$thinning,
                control = trendfilter.control.list(max_iter = obj$max_iter,
                                                   obj_tol = obj$obj_tol
                                                   )
                )
  
  obj$df <- out$df
  obj$df.min <- out$df[obj$i.min]
  obj$df.1se <- out$df[obj$i.1se]
  
  obj$tf.estimate <- glmgen:::predict.trendfilter(out,
                                                  lambda = obj$lambda.min,
                                                  x.new = obj$x.eval
                                                  ) %>% as.numeric
  
  obj$fitted.values <- glmgen:::predict.trendfilter(out,
                                                    lambda = obj$lambda.min,
                                                    x.new = obj$x
                                                    ) %>% as.numeric
  obj$residuals <- obj$y - obj$fitted.values

  obj <- obj[c("x.eval","tf.estimate","validation.method","V","validation.error.type","lambda",
               "error","se.error","lambda.min","lambda.1se","df","df.min","df.1se",
               "i.min","i.1se","x","y","weights","fitted.values","residuals","k",
               "thinning","max_iter","obj_tol")]
  
  return(obj)
}


trendfilter.validate <- function(validation.index,
                                 data.folded,
                                 obj
                                 )
  {
  
  data.train <- data.folded[-validation.index] %>% bind_rows
  data.validate <- data.folded[[validation.index]]
  
  out <- trendfilter(x = data.train$x,
                     y = data.train$y,
                     weights = data.train$weights,
                     k = obj$k,
                     lambda = obj$lambda, 
                     thinning = obj$thinning,
                     control = trendfilter.control.list(max_iter = obj$max_iter,
                                                        obj_tol = obj$obj_tol
                                                        )
                     )
  
  tf.validate.preds <- glmgen:::predict.trendfilter(out,
                                                    lambda = obj$lambda,
                                                    x.new = data.validate$x
                                                    ) %>%
    suppressWarnings
  
  if ( obj$validation.error.type == "MSE" ){
    validation.error.mat <- (tf.validate.preds - data.validate$y) ^ 2 * data.validate$weights / sum(data.validate$weights)
  }
  if ( obj$validation.error.type == "MAD" ){
    validation.error.mat <- abs(tf.validate.preds - data.validate$y) * sqrt(data.validate$weights) / sum(sqrt(data.validate$weights))
  }
  
  validation.error.sum <- colMeans(validation.error.mat) %>% as.numeric
  
  return(validation.error.sum)
}
