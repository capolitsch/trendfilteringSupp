#' Optimize the trend filtering hyperparameter (by V-fold cross validation)
#'
#' @description \loadmathjax{} \code{cv.trendfilter} performs V-fold cross 
#' validation to estimate the random-input squared error of a trend filtering 
#' estimator over a grid of hyperparameter values and returns the full error 
#' curve and the optimized trend filtering estimate within a larger list of 
#' helpful ancillary information.
#' @param x The vector of observed values of the input variable (a.k.a. the 
#' predictor, covariate, explanatory variable, regressor, independent variable, 
#' control variable, etc.)
#' @param y The vector of observed values of the output variable (a.k.a. the
#' response, target, outcome, regressand, dependent variable, etc.).
#' @param weights A vector of weights for the observed outputs. These are 
#' defined as \code{weights = 1 / sigma^2}, where \code{sigma} is a vector of 
#' standard errors of the uncertainty in the output measurements. \code{weights} 
#' should either have length equal to 1 (corresponding to observations with a 
#' constant (scalar) variance of \code{sigma = 1/sqrt(weights)}) or length equal 
#' to \code{length(y)} (i.e. heteroskedastic outputs).
#' @param k The degree of the trend filtering estimator. Defaults to 
#' \code{k=2} (quadratic trend filtering). Must be one of \code{k = 0,1,2,3},
#' although \code{k=3} is discouraged due to algorithmic instability (and is
#' visually indistinguishable from \code{k=2} anyway).
#' @param nlambda The number of trend filtering hyperparameter values 
#' to run the grid search over. If \code{lambda = NULL}, a grid of length 
#' \code{nlambda} is constructed internally.
#' @param lambda Overrides \code{nlambda} if passed. A user-supplied vector of 
#' trend filtering hyperparameter values to run the grid search over. Usually
#' good to let them be equally-spaced in log-space. Descending order is 
#' encouraged to avoid ambiguity, since they will be resorted internally for
#' algorithmic benefits. Regardless, the sorted vector is returned in the 
#' function output. Using this argument is discouraged unless you know what you 
#' are doing.
#' @param V The number of folds the data are split into for the V-fold cross
#' validation. Defaults to \code{V=5} (recommended).
#' @param validation.error.type Type of error to optimize during cross
#' validation. One of \code{c("WMAE","WMSE","MAE","MSE")}, i.e. mean-absolute 
#' deviations error, mean-squared error, and their weighted counterparts. 
#' If \code{weights = NULL}, then the weighted and 
#' unweighted counterparts are equivalent. In short, weighting helps combat
#' heteroskedasticity and absolute error decreases sensitivity to outliers.
#' Defaults to \code{"WMAE"}.
#' @param n.eval The length of the equally-spaced input grid to evaluate the 
#' optimized trend filtering estimate on.
#' @param x.eval Overrides \code{n.eval} if passed. A user-supplied grid of 
#' inputs to evaluate the optimized trend filtering estimate on. 
#' @param lambda.choice One of \code{c("lambda.min","lambda.1se")}. The choice
#' of hyperparameter that is used for optimized trend filtering estimate. \cr \cr
#' \code{lambda.min}: the hyperparameter value that minimizes the cross
#' validation error curve. \cr \cr
#' \code{lambda.1se}: the largest hyperparameter value with a cross
#' validation error within 1 standard error of the minimum cross validation 
#' error. This choice therefore favors simpler (i.e. smoother) trend filtering 
#' estimates. The motivation here is essentially Occam's razor: the two models
#' yield results that are quantitatively very close, so we favor the simpler
#' model.
#' @param thinning logical. If \code{TRUE}, then the data are preprocessed so 
#' that a smaller, better conditioned data set is used for fitting. When set to
#' \code{NULL}, the default, function will auto detect whether thinning should 
#' be applied (i.e., cases in which the numerical fitting algorithm will 
#' struggle to converge).
#' @param max_iter Maximum iterations allowed for the trend filtering 
#' convex optimization 
#' (\href{http://www.stat.cmu.edu/~ryantibs/papers/fasttf.pdf}{Ramdas and 
#' Tibshirani 2016}). Defaults to \code{max_iter = 600L}. Increase this if the 
#' trend filtering estimate does not appear to have fully converged to a 
#' reasonable estimate of the signal.
#' @param obj_tol The tolerance used in the convex optimization stopping 
#' criterion; when the relative change in the objective function is less than 
#' this value, the algorithm terminates. Decrease this if the trend 
#' filtering estimate does not appear to have fully converged to a reasonable 
#' estimate of the signal.
#' @param mc.cores Multi-core computing (for speedups): The number of cores to
#' utilize. Defaults to \code{min(detected.cores, V)}.
#' @return An object of class 'cv.trendfilter'. This is a list with the 
#' following elements:
#' \item{x.eval}{The grid of inputs the optimized trend filtering estimate was 
#' evaluated on.}
#' \item{tf.estimate}{The optimized trend filtering estimate of the signal, 
#' evaluated on \code{x.eval}.}
#' \item{validation.method}{\code{paste0(V,"-fold CV")}}
#' \item{V}{The number of folds the data are split into for the V-fold cross
#' validation.}
#' \item{validation.error.type}{Type of error that validation was performed on. 
#' One of \code{c("WMAE","WMSE","MAE","MSE")}.}
#' \item{lambda}{Vector of hyperparameter values tested during validation. This
#' vector will always be returned in descending order, regardless of the 
#' ordering provided by the user. The indices \code{i.min} and \code{i.1se}
#' correspond to this descending ordering.}
#' \item{lambda.min}{Hyperparameter value that minimizes the SURE error curve.}
#' \item{lambda.1se}{The largest hyperparameter value that is still within one
#' standard error of the minimum hyperparameter's cross validation error.}
#' \item{lambda.choice}{One of \code{c("lambda.min","lambda.1se")}. The choice
#' of hyperparameter that is used for optimized trend filtering estimate.}
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
#' \item{error}{Vector of cross validation errors for the given hyperparameter 
#' values.}
#' \item{se.error}{The standard errors of the cross validation errors.
#' These are particularly useful for implementing the ``1-standard-error rule''. 
#' The 1-SE rule favors a smoother trend filtering estimate by, instead of 
#' using the hyperparameter that minimizes the CV error, instead uses the 
#' largest hyperparameter that has a CV error within 1 standard error of the
#' smallest CV error.}
#' \item{x}{The vector of the observed inputs.}
#' \item{y}{The vector of the observed outputs.}
#' \item{weights}{A vector of weights for the observed outputs. These are
#' defined as \code{weights = 1/sigma^2}, where \code{sigma} is a vector of 
#' standard errors of the uncertainty in the measured outputs.}
#' \item{fitted.values}{The trend filtering estimate of the signal, evaluated at
#' the observed inputs \code{x}.}
#' \item{residuals}{\code{residuals = y - fitted.values}.}
#' \item{k}{The degree of the trend filtering estimator.}
#' \item{thinning}{logical. If \code{TRUE}, then the data are preprocessed so 
#' that a smaller, better conditioned data set is used for fitting. When set to
#' \code{NULL}, the default, function will auto detect whether thinning should 
#' be applied (i.e., cases in which the numerical fitting algorithm will 
#' struggle to converge).}
#' \item{max_iter}{Maximum iterations allowed for the trend filtering 
#' convex optimization 
#' (\href{http://www.stat.cmu.edu/~ryantibs/papers/fasttf.pdf}{Ramdas and 
#' Tibshirani 2016}). Increase this if the trend filtering estimate does not 
#' appear to have fully converged to a reasonable estimate of the signal.}
#' \item{obj_tol}{The tolerance used in the convex optimization stopping 
#' criterion; when the relative change in the objective function is less than 
#' this value, the algorithm terminates. Decrease this if the trend filtering 
#' estimate does not appear to have fully converged to a reasonable estimate of 
#' the signal.}
#' @details This will be a very detailed description... \cr \cr
#' \mjeqn{WMAE(\lambda) = \frac{1}{n}\sum_{i=1}^{n} |Y_i - \widehat{f}(x_i; \lambda)|\frac{\sqrt{w_i}}{\sum_j\sqrt{w_j}}}{ascii} \cr 
#' \mjeqn{WMSE(\lambda) = \frac{1}{n}\sum_{i=1}^{n} |Y_i - \widehat{f}(x_i; \lambda)|^2\frac{w_i}{\sum_jw_j}}{ascii} \cr 
#' \mjeqn{MAE(\lambda) = \frac{1}{n}\sum_{i=1}^{n} |Y_i - \widehat{f}(x_i; \lambda)|}{ascii} \cr 
#' \mjeqn{MSE(\lambda) = \frac{1}{n}\sum_{i=1}^{n} |Y_i - \widehat{f}(x_i; \lambda)|^2}{ascii} \cr \cr 
#' where \mjeqn{\widehat{f}(x_i; \lambda)}{ascii} is the trend filtering 
#' estimate with hyperparameter \eqn{\lambda}, evaluated at 
#' \mjeqn{x_i}{ascii}. 
#' @export cv.trendfilter
#' @author Collin A. Politsch, \email{collinpolitsch@@gmail.com}
#' @seealso \code{\link{SURE.trendfilter}}, \code{\link{bootstrap.trendfilter}}
#' @references 
#' \strong{Cross validation}
#' \enumerate{
#' \item Hastie, Tibshirani, and Friedman (2009). The Elements of Statistical 
#' Learning: Data Mining, Inference, and Prediction. 2nd edition. Springer 
#' Series in Statistics. 
#' \href{https://web.stanford.edu/~hastie/ElemStatLearn/printings/ESLII_print12_toc.pdf}{
#' [Online print #12]}. (See Sections 7.10 and 7.12) \cr
#' \item James, Witten, Hastie, and Tibshirani (2013). An Introduction to 
#' Statistical Learning : with Applications in R. Springer.
#' \href{https://www.statlearning.com/}{[Most recent online print]} (See 
#' Section 5.1). \emph{Less technical than the above reference.}\cr
#' \item Tibshirani (2013). Model selection and validation 2: Model
#' assessment, more cross-validation. \emph{36-462: Data Mining course notes} 
#' (Carnegie Mellon).
#' \href{https://www.stat.cmu.edu/~ryantibs/datamining/lectures/19-val2.pdf}{[Link]}
#' }
#' \strong{Trend filtering optimization algorithm}
#' \enumerate{
#' \item{Ramdas and Tibshirani (2016). Fast and Flexible ADMM Algorithms 
#' for Trend Filtering. \emph{Journal of Computational and Graphical 
#' Statistics}, 25(3), p. 839-858.
#' \href{https://amstat.tandfonline.com/doi/abs/10.1080/10618600.2015.1054033#.XfJpNpNKju0}{[Link]}} \cr
#' \item{Arnold, Sadhanala, and Tibshirani (2014). Fast algorithms for 
#' generalized lasso problems. R package \emph{glmgen}. Version 0.0.3. 
#' \href{https://github.com/glmgen/glmgen}{[Link]}} \cr
#' (Software implementation of Ramdas and Tibshirani algorithm) \cr
#' }
#' @examples 
#' #######################################################################
#' ###  Phase-folded light curve of an eclipsing binary star system   ####
#' #######################################################################
#' 
#' # A binary star system is a pair of closely-separated stars that move
#' # in an orbit around a common center of mass. When a binary star system 
#' # is oriented in such a way that the stars eclipse one another from our 
#' # vantage point on Earth, we call it an 'eclipsing binary (EB) star 
#' # system'. From our perspective, the total brightness of an EB dips 
#' # periodically over time due to the stars eclipsing one another. And 
#' # the shape of the brightness curve is consistent within each period
#' # of the orbit. In order to learn about the physics of these EBs,
#' # astronomers 'phase-fold' the brightness curve so that all the orbital 
#' # periods are stacked on top of one another in a plot of the EB's phase 
#' # vs. its apparent brightness, and then find a 'best-fitting' model
#' # for the phase-folded curve. Here, we use trend filtering to fit an
#' # optimal phase-folded model for an EB.
#' 
#' data(eclipsing_binary)
#' 
#' # head(data)
#' #
#' # |      phase|      flux|  std.err|
#' # |----------:|---------:|--------:|
#' # | -0.4986308| 0.9384845| 0.010160|
#' # | -0.4978067| 0.9295757| 0.010162|
#' # | -0.4957892| 0.9438493| 0.010162|
#' 
#' cv.out <- cv.trendfilter(x = data$phase, 
#'                          y = data$flux, 
#'                          weights = 1 / data$std.err ^ 2,
#'                          validation.error.type = "MAE",
#'                          thinning = TRUE, 
#'                          obj_tol = 1e-14, 
#'                          max_iter = 1e4)
#' 
#' # Plot the results
#' 
#' par(mfrow = c(2,1), mar = c(5,4,2.5,1) + 0.1)
#' plot(log(cv.out$lambda), cv.out$error, main = "CV error curve", 
#'      xlab = "log(lambda)", ylab = "CV error")
#' segments(x0 = log(cv.out$lambda), x1 = log(cv.out$lambda), 
#'          y0 = cv.out$error - cv.out$se.error, 
#'          y1 = cv.out$error + cv.out$se.error)
#' abline(v = log(cv.out$lambda.min), lty = 2, col = "blue3")
#' text(x = log(cv.out$lambda.min), y = par("usr")[4], 
#'      labels = "optimal hyperparameter", pos = 1, col = "blue3")
#' plot(data$phase, data$flux, cex = 0.15, xlab = "Phase", ylab = "Flux",
#'      main = "Eclipsing binary phase-folded light curve")
#' segments(x0 = data$phase, x1 = data$phase, 
#'          y0 = data$flux - data$std.err, y1 = data$flux + data$std.err, 
#'          lwd = 0.25)
#' lines(cv.out$x.eval, cv.out$tf.estimate, col = "orange", lwd = 2.5)


#' @importFrom dplyr arrange case_when group_split bind_rows
#' @importFrom magrittr %$% %>%
#' @importFrom tidyr drop_na
#' @importFrom parallel mclapply detectCores
#' @importFrom matrixStats rowSds
cv.trendfilter <- function(x, 
                           y, 
                           weights = NULL, 
                           k = 2L, 
                           nlambda = 250L,
                           lambda = NULL, 
                           validation.error.type = c("WMAE","WMSE","MAE","MSE"),
                           V = 5L,
                           lambda.choice = c("lambda.min","lambda.1se"),
                           n.eval = 1500L,
                           x.eval = NULL,
                           thinning = NULL,
                           max_iter = 600L, 
                           obj_tol = 1e-10,
                           mc.cores = detectCores()
                           )
  {

  if ( missing(x) || is.null(x) ) stop("x must be passed.")
  if ( missing(y) || is.null(y) ) stop("y must be passed.")
  if ( length(x) != length(y) ) stop("x and y must have the same length.")
  if ( !is.null(weights) ){
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
  if ( is.null(x.eval) ){
    if ( n.eval != round(n.eval) ) stop("n.eval must be a positive integer.")
  }else{
    if ( any(x.eval < min(x) | x.eval > max(x)) ) stop("x.eval should all be in range(x).")
  }
  if ( mc.cores > detectCores() ){
    warning(paste0("Your machine only has ", detectCores(), " cores. Adjusting mc.cores accordingly."))
    mc.cores <- detectCores()
  }
  
  mc.cores <- min(mc.cores, V)
  lambda.choice <- match.arg(lambda.choice)
  validation.error.type <- match.arg(validation.error.type)
  
  if ( length(weights) != length(y) ){
    if ( length(weights) == 0 ){
      weights <- rep_len(1, length(y))
    }
    if ( length(weights) == 1 ){
      weights <- rep_len(weights, length(y))
    }
  }
  
  data <- tibble(x, y, weights) %>% 
    arrange(x) %>% 
    filter( weights != 0 ) %>%
    drop_na 
  
  if ( is.null(lambda) ){
    lambda <- seq(16, -10, length = nlambda) %>% exp 
  }else{
    lambda <- sort(lambda, decreasing = TRUE)
  }
  
  x.scale <- mean(diff(data$x))
  y.scale <- mean(abs(data$y)) / 10
  x <- data$x / x.scale
  y <- data$y / y.scale
  if ( !all(weights == 1) ){
    weights <- y.scale ^ 2 * data$weights
  }
  
  data.folded <- tibble(x, y, weights) %>% 
    group_split( sample( rep_len(1:V, nrow(data)) ), .keep = FALSE )
  
  if ( is.null(x.eval) ){
    x.eval <- seq(min(x), max(x), length = n.eval)
  }else{
    x.eval <- sort(x.eval) / x.scale
  }
  
  obj <- structure(list(x.eval = x.eval,
                        validation.method = paste0(V,"-fold CV"),
                        V = as.integer(V),
                        validation.error.type = validation.error.type,
                        lambda = lambda, 
                        lambda.choice = lambda.choice,
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
  
  rm(x,y,weights,lambda,k,thinning,max_iter,obj_tol,data)

  cv.out <- matrix(unlist(mclapply(1:(obj$V), 
                                   FUN = trendfilter.validate, 
                                   data.folded = data.folded, 
                                   obj = obj, 
                                   mc.cores = mc.cores
                                   )
                          ), 
                   ncol = obj$V
                   )
  
  obj <- c(obj, list(error = rowMeans(cv.out),
                     se.error = rowSds(cv.out) / sqrt(obj$V)
                     )
           )
 
  obj$i.min <- which.min(obj$error)
  obj$lambda.min <- obj$lambda[obj$i.min]
  obj$i.1se <- obj %$% which(error <= error[i.min] + se.error[i.min]) %>% min
  obj$lambda.1se <- obj$lambda[obj$i.1se]
  
  optimization.controls <- glmgen::trendfilter.control.list(max_iter = obj$max_iter,
                                                            obj_tol = obj$obj_tol
                                                            )

  out <- obj %$%
    glmgen::trendfilter(x = x,
                        y = y,
                        weights = weights,
                        lambda = lambda,
                        k = k,
                        thinning = thinning,
                        control = optimization.controls
                        )
  
  obj$df <- out$df
  obj$df.min <- out$df[obj$i.min]
  obj$df.1se <- out$df[obj$i.1se]
  
  lambda.pred <- case_when(
    obj$lambda.choice == "lambda.min" ~ obj$lambda.min,
    obj$lambda.choice == "lambda.1se" ~ obj$lambda.1se
  )
  
  obj$tf.estimate <- glmgen:::predict.trendfilter(out,
                                                  lambda = lambda.pred,
                                                  x.new = obj$x.eval
                                                  ) %>% as.numeric
  
  obj$fitted.values <- glmgen:::predict.trendfilter(out,
                                                    lambda = lambda.pred,
                                                    x.new = obj$x
                                                    ) %>% as.numeric

  obj$x.eval <- obj$x.eval * x.scale
  obj$tf.estimate <- obj$tf.estimate * y.scale
  obj$x <- obj$x * x.scale
  obj$y <- obj$y * y.scale
  obj$weights <- obj$weights / y.scale ^ 2
  obj$fitted.values <- obj$fitted.values * y.scale
  obj$residuals <- (obj$y - obj$fitted.values) * y.scale

  obj <- obj[c("x.eval","tf.estimate","validation.method","V","validation.error.type",
               "lambda","lambda.min","lambda.1se","lambda.choice","error","se.error",
               "df","df.min","df.1se","i.min","i.1se","x","y","weights","fitted.values",
               "residuals","k","thinning","max_iter","obj_tol")]
  
  return(obj)
}


trendfilter.validate <- function(validation.index,
                                 data.folded,
                                 obj
                                 )
  {
  
  data.train <- data.folded[-validation.index] %>% bind_rows
  data.validate <- data.folded[[validation.index]]
  
  optimization.controls <- glmgen::trendfilter.control.list(max_iter = obj$max_iter,
                                                            obj_tol = obj$obj_tol
                                                            )
  
  out <- glmgen::trendfilter(x = data.train$x,
                             y = data.train$y,
                             weights = data.train$weights,
                             k = obj$k,
                             lambda = obj$lambda,
                             thinning = obj$thinning,
                             control = optimization.controls
                             )
  
  tf.validate.preds <- glmgen:::predict.trendfilter(out,
                                                    lambda = obj$lambda,
                                                    x.new = data.validate$x
                                                    ) %>%
    suppressWarnings
  
  if ( obj$validation.error.type == "MSE" ){
    validation.error.mat <- obj$y.scale ^ 2 * (tf.validate.preds - data.validate$y) ^ 2
  }
  if ( obj$validation.error.type == "MAE" ){
    validation.error.mat <- obj$y.scale * abs(tf.validate.preds - data.validate$y)
  }
  if ( obj$validation.error.type == "WMSE" ){
    validation.error.mat <- obj$y.scale ^ 2 * (tf.validate.preds - data.validate$y) ^ 2 * 
    (data.validate$weights) / sum((data.validate$weights))
  }
  if ( obj$validation.error.type == "WMAE" ){
    validation.error.mat <- obj$y.scale * abs(tf.validate.preds - data.validate$y) * 
      sqrt(data.validate$weights) / sum(sqrt(data.validate$weights))
  }
  
  validation.error.sum <- colMeans(validation.error.mat) %>% as.numeric
  
  return(validation.error.sum)
}
