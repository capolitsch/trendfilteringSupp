#' Optimize the trend filtering hyperparameter (by V-fold cross validation)
#'
#' @description \loadmathjax{} \code{cv.trendfilter} performs V-fold cross validation to
#' estimate the random-input squared error of a trend filtering estimator over 
#' a grid of hyperparameter values and returns the optimized estimator.
#' @param x The vector of the observed inputs. 
#' \mjeqn{x+y}{ascii}
#' @param y The vector of the observed outputs.
#' @param weights A vector of weights for the observed outputs. These are defined as 
#' \code{weights = 1 / sigma^2}, where \code{sigma} is a vector of standard 
#' errors of the uncertainty in the output measurements. \code{weights} should 
#' either have length equal to 1 (corresponding to observations with a constant 
#' (scalar) variance of \code{sigma = 1/sqrt(weights)}) or length equal to 
#' \code{length(y)} (i.e. heteroskedastic outputs).
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
#' @param V The number of folds the data are split into for the V-fold cross
#' validation. Defaults to \code{V=10} (recommended).
#' \code{V=length(x)} is equivalent to leave-one-out cross validation.
#' @param validation.error.type Type of error to optimize during cross
#' validation. One of c("WMAD","MAD","WMSE","MSE"), i.e. mean-absolute 
#' deviations error, mean-squared error, and their weighted counterparts. If 
#' \code{weights = NULL}, then each weighted and non-weighted pair are 
#' equivalent. Defaults to \code{"WMAD"}. That is,
#' \deqn{WMAD(\lambda) = 1/n \sum |y - tf.estimate_i| * \sigma}
#' @param n.eval (integer) The length of the equally-spaced input grid to 
#' evaluate the optimized trend filtering estimate on.
#' @param x.eval Overrides \code{n.eval} if passed. A user-supplied grid of 
#' inputs to evaluate the optimized trend filtering estimate on. 
#' @param thinning (logical) If \code{TRUE}, then the data are preprocessed so 
#' that a smaller, better conditioned data set is used for fitting. When set to
#' \code{NULL}, the default, function will auto detect whether thinning should 
#' be applied (i.e., cases in which the numerical fitting algorithm will 
#' struggle to converge).
#' @param max_iter Maximum iterations allowed for the trend filtering 
#' convex optimization 
#' [\href{http://www.stat.cmu.edu/~ryantibs/papers/fasttf.pdf}{Ramdas & Tibshirani (2015)}]. 
#' Defaults to \code{max_iter = 250}. Consider increasing this if the trend 
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
#' \item{lambda}{Vector of hyperparameter values tested during validation. This
#' vector will always be returned in descending order, regardless of the 
#' ordering provided by the user. The indices \code{i.min} and \code{i.1se}
#' correspond to this descending ordering.}
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
#' \item \href{https://web.stanford.edu/~hastie/ElemStatLearn/printings/ESLII_print12_toc.pdf}{
#' Hastie et al. (2017). The Elements of Statistical Learning}
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
#' # Run the cross validation for a quadratic trend filtering estimator, i.e. 
#' # k = 2 (default)
#' 
#' cv.obj <- cv.trendfilter(x = log10.wavelength,
#'                          y = flux,
#'                          weights = weights)
#' 
#'                                           
#' # Extract the CV error curve and optimized trend filtering estimate from 
#' # the output
#' 
#' log.lambda <- log(cv.obj$lambda)
#' error <- cv.obj$error
#' log.lambda.min <- log(cv.obj$lambda.min)
#' 
#' log10.wavelength.eval <- cv.obj$x.eval
#' tf.estimate <- cv.obj$tf.estimate
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
#' plot(x = log.lambda, y = error, main = "CV error curve", 
#'      xlab = "log(lambda)", ylab = "CV error")
#' abline(v = log.lambda.min, lty = 2, col = "blue3")
#' text(x = log.lambda.min, y = par("usr")[4], 
#'      labels = "optimal hyperparameter", pos = 1, col = "blue3")
#' 
#' plot(x = wavelength, y = flux, type = "l", main = "Quasar Lyman-alpha forest", 
#'      xlab = "Observed wavelength (Angstroms)", ylab = "Flux")
#' lines(wavelength.eval, tf.estimate, col = "orange", lwd = 2.5)
#' legend(x = "topleft", lwd = c(1,2), lty = 1, col = c("black","orange"), 
#'        legend = c("Noisy quasar Lyman-alpha forest", "Trend filtering estimate"))

#' @import magrittr
#' @import glmgen
#' @importFrom parallel mclapply
#' @importFrom matrixStats rowSds
cv.trendfilter <- function(x, 
                           y, 
                           weights = NULL, 
                           k = 2L, 
                           nlambda = 250L,
                           lambda = NULL, 
                           V = 10L,
                           validation.error.type = c("MAD","WMAD","MSE","WMSE"),
                           n.eval = 1500L,
                           x.eval = NULL,
                           thinning = NULL,
                           max_iter = 500L, 
                           obj_tol = 1e-10,
                           mc.cores = max(c(parallel::detectCores() - 2), 1)
                           )
  {

  if ( missing(x) || is.null(x) ) stop("x must be passed.")
  if ( missing(y) || is.null(y) ) stop("y must be passed.")
  if ( length(x) != length(y) ) stop("x and y must have the same length.")
  if ( !is.null(weights) ){
    if ( any(weights == 0L) ) stop("cannot pass zero weights.")
    if ( !(length(weights) %in% c(1,length(y))) ) stop("weights must either be have length 1 or length(y), or be NULL.")
  }
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
  
  validation.error.type <- match.arg(validation.error.type)
  mc.cores <- ifelse(mc.cores > V, V, mc.cores)
  
  if ( is.null(lambda) ){
    lambda <- seq(10, -10, length = nlambda) %>% exp 
  }else{
    lambda <- sort(lambda, decreasing = TRUE)
  }

  if ( length(weights) != length(y) ){
    weights <- case_when(
      length(weights) == 1 ~ rep_len(weights, length(y)),
      length(weights) == 0 ~ rep_len(1, length(y))
    )
  }

  data <- tibble(x, y, weights) %>% 
    arrange(x) %>% 
    drop_na
  
  x.scale <- mean(diff(data$x))
  y.scale <- mean(data$y) / 10
  x <- data$x / x.scale
  y <- data$y / y.scale
  weights <- y.scale ^ 2 * data$weights
  
  data.folded <- tibble(x, y, weights) %>% 
    group_split( sample( rep_len(1:V, nrow(data)) ), .keep = FALSE )
  
  if ( is.null(x.eval) ){
    x.eval <- seq(min(x), max(x), length = n.eval)
  }else{
    x.eval <- sort(x.eval) / x.scale
  }
  
  obj <- structure(list(x.eval = x.eval,
                        validation.method = "cv",
                        V = as.integer(V),
                        validation.error.type = validation.error.type,
                        lambda = lambda, 
                        x = x,
                        y = y,
                        y.scale = y.scale,
                        weights = weights,
                        k = as.integer(k),
                        thinning = thinning,
                        max_iter = as.integer(max_iter),
                        obj_tol = obj_tol
                        ),
                   class = "cv.trendfilter"
                   )
  
  rm(x,y,weights,lambda,k,thinning,max_iter,obj_tol,V,data)

  cv.out <- mclapply(1:obj$V, 
                     FUN = trendfilter.validate, 
                     data.folded = data.folded, 
                     obj = obj, 
                     mc.cores = mc.cores
                     ) %>%
    unlist %>%
    matrix(ncol = obj$V)
  
  obj <- c(obj, list(error = rowMeans(cv.out),
                     se.error = matrixStats::rowSds(cv.out) / sqrt(obj$V)
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
                lambda = lambda,
                k = k, 
                thinning = thinning,
                control = trendfilter.control.list(max_iter = max_iter,
                                                   obj_tol = obj_tol
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
  obj$x.eval <- obj$x.eval * x.scale
  obj$tf.estimate <- obj$tf.estimate * y.scale
  obj$x <- obj$x * x.scale
  obj$y <- obj$y * y.scale
  obj$weights <- obj$weights * y.scale ^ 2
  obj$fitted.values <- obj$fitted.values * y.scale
  obj$residuals <- (obj$y - obj$fitted.values) * y.scale

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
    validation.error.mat <- obj$y.scale ^ 2 * (tf.validate.preds - data.validate$y) ^ 2
  }
  if ( obj$validation.error.type == "MAD" ){
    validation.error.mat <- obj$y.scale * abs(tf.validate.preds - data.validate$y)
  }
  if ( obj$validation.error.type == "WMSE" ){
    validation.error.mat <- obj$y.scale ^ 2 * (tf.validate.preds - data.validate$y) ^ 2 * 
    (data.validate$weights / obj$y.scale ^ 2) / sum((data.validate$weights) / obj$y.scale ^ 2)
  }
  if ( obj$validation.error.type == "WMAD" ){
    validation.error.mat <- obj$y.scale * abs(tf.validate.preds - data.validate$y) * 
      sqrt(data.validate$weights / obj$y.scale ^ 2) / sum(sqrt(data.validate$weights / obj$y.scale ^ 2))
  }
  
  validation.error.sum <- colMeans(validation.error.mat) %>% as.numeric
  
  return(validation.error.sum)
}
