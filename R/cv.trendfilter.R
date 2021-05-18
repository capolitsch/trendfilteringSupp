#' Optimize the trend filtering hyperparameter (by V-fold cross validation)
#'
#' @description \loadmathjax{} \code{cv.trendfilter} performs V-fold cross 
#' validation to estimate the random-input squared error of a trend filtering 
#' estimator on a grid of values for the hyperparameter \code{gamma}, and 
#' returns the full error curve and the optimized trend filtering estimate 
#' within a larger list with useful ancillary information.
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
#' to \code{length(y)} (i.e. general heteroskedastic noise).
#' @param k The degree of the trend filtering estimator. Defaults to 
#' \code{k=2} (quadratic trend filtering). Must be one of \code{k = 0,1,2,3},
#' although \code{k=3} is discouraged due to algorithmic instability (and is
#' visually indistinguishable from \code{k=2} anyway).
#' @param ngammas Integer. The number of trend filtering hyperparameter values 
#' to run the grid search over.
#' @param gammas Overrides \code{ngammas} if passed. A user-supplied vector of 
#' trend filtering hyperparameter values to run the grid search over. It is
#' advisable to let the vector be equally-spaced in log-space and provided in 
#' descending order. The function output will contain the sorted hyperparameter
#' vector regardless of the input ordering, and all related output objects 
#' (e.g. the \code{errors} vector) will correspond to the sorted ordering. 
#' Unless you know what you are doing, it is best to leave this \code{NULL}.
#' validation. 
#' @param V The number of folds the data are split into for the V-fold cross
#' validation. Defaults to \code{V=5} (recommended).
#' @param validation.error.type Type of error to optimize during cross
#' validation. One of \code{c("WMAE","WMSE","MAE","MSE")}, i.e. mean-absolute 
#' deviations error, mean-squared error, and their weighted counterparts. 
#' If \code{weights = NULL}, then the weighted and 
#' unweighted counterparts are equivalent. In short, weighting helps combat
#' heteroskedasticity and absolute error decreases sensitivity to outliers.
#' Defaults to \code{"WMAE"}.
#' @param nx.eval The length of the equally-spaced input grid to evaluate the 
#' evaluate the optimized trend filtering estimate on.
#' @param x.eval Overrides \code{nx.eval} if passed. A user-supplied grid of 
#' inputs to evaluate the optimized trend filtering estimate on. 
#' @param gamma.choice One of \code{c("gamma.min","gamma.1se")}. The choice
#' of hyperparameter that is used for optimized trend filtering estimate. 
#' \cr \cr
#' \code{gamma.min}: the hyperparameter value that minimizes the cross
#' validation error curve. \cr \cr
#' \code{gamma.1se}: the largest hyperparameter value with a cross
#' validation error within 1 standard error of the minimum cross validation 
#' error. This choice therefore favors simpler (i.e. smoother) trend filtering 
#' estimates. The motivation here is essentially Occam's razor: the two models
#' yield results that are quantitatively very close, so we favor the simpler
#' model.
#' @param thinning logical. If \code{TRUE}, then the data are preprocessed so 
#' that a smaller, better conditioned data set is used for fitting. When left
#' \code{NULL}, the default, the optimization will automatically detect whether 
#' thinning should be applied (i.e., cases in which the numerical fitting 
#' algorithm will struggle to converge). This preprocessing procedure is 
#' controlled by the \code{x_tol} argument of 
#' \code{\link[glmgen]{trendfilter.control.list}}.
#' @param optimization.params a named list of parameters produced by the
#' \pkg{glmgen} function \code{\link[glmgen]{trendfilter.control.list}} that
#' contains all parameter choices (user-supplied or defaults) to be passed to 
#' the trend filtering ADMM algorithm
#' (\href{http://www.stat.cmu.edu/~ryantibs/papers/fasttf.pdf}{Ramdas and
#' Tibshirani 2016}). See the linked function documentation for full details. 
#' No technical understanding of the ADMM algorithm is needed and the default
#' parameter choices will almost always suffice. However, the following three 
#' parameters may require some adjustments to ensure that your trend filtering
#' estimate has sufficiently converged:
#' \enumerate{ 
#' \item{\code{max_iter}}: Maximum iterations allowed for the trend filtering 
#' convex optimization. Defaults to \code{max_iter = 600L}. Increase this if 
#' the trend filtering estimate does not appear to have fully converged to a 
#' reasonable estimate of the signal.
#' \item{\code{obj_tol}}: The tolerance used in the convex optimization stopping 
#' criterion; when the relative change in the objective function is less than 
#' this value, the algorithm terminates. Decrease this if the trend filtering 
#' estimate does not appear to have fully converged to a reasonable estimate of 
#' the signal.
#' \item{x_tol}: defines uniqueness or sameness of \code{x}'s. If we make bins 
#' of size \code{x_tol} and find at least two \code{x}'s which fall into the 
#' same bin, then we thin the data.
#' }
#' @param mc.cores Multi-core computing (for speedups): The number of cores to
#' utilize. Defaults to the number of cores detected.
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
#' \item{gammas}{Vector of hyperparameter values tested during validation. This
#' vector will always be returned in descending order, regardless of the 
#' ordering provided by the user. The indices \code{i.min} and \code{i.1se}
#' correspond to this descending ordering.}
#' \item{gamma.min}{Hyperparameter value that minimizes the SURE error curve.}
#' \item{gamma.1se}{The largest hyperparameter value that is still within one
#' standard error of the minimum hyperparameter's cross validation error.}
#' \item{gamma.choice}{One of \code{c("gamma.min","gamma.1se")}. The choice
#' of hyperparameter that is used for optimized trend filtering estimate.}
#' \item{edfs}{Vector of effective degrees of freedom for trend filtering
#' estimators fit during validation.}
#' \item{edf.min}{The effective degrees of freedom of the optimally-tuned trend 
#' filtering estimator.}
#' \item{edf.1se}{The effective degrees of freedom of the 1-stand-error rule
#' trend filtering estimator.}
#' \item{i.min}{The index of \code{gammas} that minimizes the cross validation 
#' error.}
#' \item{i.1se}{The index of \code{gammas} that gives the largest hyperparameter
#' value that has a cross validation error within 1 standard error of the 
#' minimum of the cross validation error curves.}
#' \item{errors}{Vector of cross validation errors for the given hyperparameter 
#' values.}
#' \item{se.errors}{The standard errors of the cross validation errors.
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
#' that a smaller, better conditioned data set is used for fitting.}
#' \item{optimization.params}{a list of parameters that control the trend
#' filtering convex optimization.}
#' \item{n.iter}{Vector of the number of iterations needed for the ADMM
#' algorithm to converge within the given tolerance, for each hyperparameter
#' value. If many of these are exactly equal to \code{max_iter}, then their
#' solutions have not converged with the tolerance specified by \code{obj_tol}.
#' In which case, it is often prudent to increase \code{max_iter}.}
#' \item{x.scale, y.scale, data.scaled}{for internal use.}
#' @details This will be a very detailed description... \cr \cr
#' \mjeqn{WMAE(\gamma) = \frac{1}{n}\sum_{i=1}^{n} |Y_i - \widehat{f}(x_i; \gamma)|\frac{\sqrt{w_i}}{\sum_j\sqrt{w_j}}}{ascii} \cr 
#' \mjeqn{WMSE(\gamma) = \frac{1}{n}\sum_{i=1}^{n} |Y_i - \widehat{f}(x_i; \gamma)|^2\frac{w_i}{\sum_jw_j}}{ascii} \cr 
#' \mjeqn{MAE(\gamma) = \frac{1}{n}\sum_{i=1}^{n} |Y_i - \widehat{f}(x_i; \gamma)|}{ascii} \cr 
#' \mjeqn{MSE(\gamma) = \frac{1}{n}\sum_{i=1}^{n} |Y_i - \widehat{f}(x_i; \gamma)|^2}{ascii} \cr \cr 
#' where \mjeqn{\widehat{f}(x_i; \gamma)}{ascii} is the trend filtering 
#' estimate with hyperparameter \eqn{\gamma}, evaluated at 
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
#' \href{https://www.stat.cmu.edu/~ryantibs/datamining/lectures/19-val2.pdf}{
#' [Link]}
#' }
#' \strong{Trend filtering optimization algorithm}
#' \enumerate{
#' \item{Ramdas and Tibshirani (2016). Fast and Flexible ADMM Algorithms 
#' for Trend Filtering. \emph{Journal of Computational and Graphical 
#' Statistics}, 25(3), p. 839-858.
#' \href{https://amstat.tandfonline.com/doi/abs/10.1080/10618600.2015.1054033#.XfJpNpNKju0}{
#' [Link]}} \cr
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
#' # head(df)
#' #
#' # |      phase|      flux|  std.err|
#' # |----------:|---------:|--------:|
#' # | -0.4986308| 0.9384845| 0.010160|
#' # | -0.4978067| 0.9295757| 0.010162|
#' # | -0.4957892| 0.9438493| 0.010162|
#' 
#' # This specific choice of values did not come a priori
#' gamma.grid <- exp( seq(7, 16, length = 150) )
#' 
#' cv.out <- cv.trendfilter(x = df$phase, 
#'                          y = df$flux, 
#'                          weights = 1 / df$std.err ^ 2,
#'                          gammas = gamma.grid,
#'                          validation.error.type = "MAE",
#'                          thinning = TRUE, 
#'                          optimization.params = trendfilter.control.list(max_iter = 5e3,
#'                                                                         obj_tol = 1e-6)
#'                          )
#' 
#' # Plot the results
#' 
#' par(mfrow = c(2,1), mar = c(5,4,2.5,1) + 0.1)
#' plot(log(cv.out$gammas), cv.out$errors, main = "CV error curve", 
#'      xlab = "log(gamma)", ylab = "CV error")
#' segments(x0 = log(cv.out$gammas), x1 = log(cv.out$gammas), 
#'          y0 = cv.out$errors - cv.out$se.errors, 
#'          y1 = cv.out$errors + cv.out$se.errors)
#' abline(v = log(cv.out$gamma.min), lty = 2, col = "blue3")
#' text(x = log(cv.out$gamma.min), y = par("usr")[4], 
#'      labels = "optimal gamma", pos = 1, col = "blue3")
#' plot(df$phase, df$flux, cex = 0.15, xlab = "Phase", ylab = "Flux",
#'      main = "Eclipsing binary phase-folded light curve")
#' segments(x0 = df$phase, x1 = df$phase, 
#'          y0 = df$flux - df$std.err, y1 = df$flux + df$std.err, 
#'          lwd = 0.25)
#' lines(cv.out$x.eval, cv.out$tf.estimate, col = "orange", lwd = 2.5)


#' @importFrom dplyr arrange case_when group_split bind_rows
#' @importFrom glmgen trendfilter.control.list
#' @importFrom parallel mclapply detectCores
#' @importFrom matrixStats rowSds
#' @importFrom magrittr %$% %>%
#' @importFrom tidyr drop_na
cv.trendfilter <- function(x, 
                           y, 
                           weights = NULL, 
                           k = 2L, 
                           V = 5L,
                           ngammas = 250L,
                           gammas = NULL, 
                           gamma.choice = c("gamma.min","gamma.1se"),
                           validation.error.type = c("WMAE","WMSE","MAE","MSE"),
                           nx.eval = 1500L,
                           x.eval = NULL,
                           thinning = NULL,
                           optimization.params = trendfilter.control.list(max_iter = 600L,
                                                                          obj_tol = 1e-10),
                           mc.cores = detectCores()
                           )
  {

  if ( missing(x) || is.null(x) ) stop("x must be passed.")
  if ( missing(y) || is.null(y) ) stop("y must be passed.")
  if ( length(x) != length(y) ) stop("x and y must have the same length.")
  if ( !is.null(weights) ){
    if ( !(length(weights) %in% c(1,length(y))) ){
      stop("weights must either be have length 1 or length(y), or be NULL.")
    }
  }
  if ( length(y) < k + 2 ){
    stop("y must have length >= k+2 for kth order trend filtering.")
  }
  if ( k < 0 || k != round(k) ){
    stop("k must be a nonnegative integer. k=2 recommended")
  }
  if ( k > 3 ){
    stop(paste0("Large k leads to generally worse conditioning.\n", 
                "k = 0,1,2 are the most stable choices."))
  }
  if ( k == 3 ){
    warning(paste0("k = 3 can have poor conditioning...\n", 
                   "k = 2 is more stable and visually indistinguishable."))
  }
  if ( is.null(gammas) && (ngammas != round(ngammas) || ngammas < 25L) ){
    stop("ngammas must be a positive integer >= 25.")
  }
  if ( !is.null(gammas) ){
    if ( min(gammas) < 0L ){
      stop("All specified gamma values must be positive.")
    }
    if ( length(gammas) < 25L ) stop("gammas must have length >= 25.")
  }
  if ( is.null(x.eval) ){
    if ( nx.eval != round(nx.eval) ) stop("nx.eval must be a positive integer.")
  }else{
    if ( any(x.eval < min(x) || x.eval > max(x)) ){
      stop("x.eval should all be in range(x).")
    }
  }
  if ( mc.cores > detectCores() ){
    warning(paste0("Your machine only has ", detectCores(), 
                   " cores. Adjusting mc.cores accordingly."))
    mc.cores <- detectCores()
  }
  if ( length(weights) == 1 ){
    weights <- rep(weights, length(y))
  }
  if ( length(weights) == 0 ){
    weights <- rep(1, length(y))
  }
  if ( is.null(gammas) ){
    gammas <- seq(16, -10, length = ngammas) %>% exp 
  }else{
    gammas <- sort(gammas, decreasing = TRUE)
  }
  
  mc.cores <- min(mc.cores, V)
  gamma.choice <- match.arg(gamma.choice)
  validation.error.type <- match.arg(validation.error.type)
  
  data <- tibble(x, y, weights) %>% 
    arrange(x) %>% 
    filter( weights != 0 ) %>%
    drop_na 
  
  rm(x,y,weights)
  
  x.scale <- median(diff(data$x))
  y.scale <- median(abs(data$y)) / 10
  optimization.params$x_tol <- optimization.params$x_tol / x.scale
  
  data.scaled <- data %>%
    mutate(x = x / x.scale,
           y = y / y.scale,
           weights = weights * y.scale ^ 2 )
  
  data.folded <- data.scaled %>% 
    group_split( sample( rep_len(1:V, nrow(data.scaled)) ), .keep = FALSE )
  
  if ( is.null(x.eval) ){
    x.eval <- seq(min(data$x), max(data$x), length = nx.eval)
  }else{
    x.eval <- sort(x.eval)
  }

  obj <- structure(list(x.eval = x.eval,
                        validation.method = paste0(V,"-fold CV"),
                        V = as.integer(V),
                        validation.error.type = validation.error.type,
                        gammas = gammas, 
                        gamma.choice = gamma.choice,
                        x = data$x,
                        y = data$y,
                        weights = data$weights,
                        k = as.integer(k),
                        thinning = thinning,
                        optimization.params = optimization.params,
                        data.scaled = data.scaled,
                        x.scale = x.scale, 
                        y.scale = y.scale
                        ),
                   class = "cv.trendfilter"
                   )
  
  rm(V,validation.error.type,gammas,ngammas,gamma.choice,k,thinning,data,nx.eval,
     optimization.params,data.scaled,x.eval,x.scale,y.scale)

  cv.out <- matrix(unlist(mclapply(1:(obj$V), 
                                   FUN = trendfilter.validate, 
                                   data.folded = data.folded, 
                                   obj = obj, 
                                   mc.cores = mc.cores
                                   )
                          ), 
                   ncol = obj$V
                   )
  
  errors <- rowMeans(cv.out)
  se.errors <- rowSds(cv.out) / sqrt(obj$V)
 
  obj$i.min <- which.min(errors)
  obj$i.1se <- which(errors <= errors[obj$i.min] + se.errors[obj$i.min]) %>% min
  obj$gamma.min <- obj$gammas[obj$i.min]
  obj$gamma.1se <- obj$gammas[obj$i.1se]
  
  if ( obj$validation.error.type %in% c("MSE","WMSE") ){
    obj$errors <- errors * obj$y.scale ^ 2
    obj$se.errors <- se.errors * obj$y.scale ^ 2
  }
  if ( obj$validation.error.type %in% c("MAE","WMAE") ){
    obj$errors <- errors * obj$y.scale
    obj$se.errors <- se.errors * obj$y.scale
  }

  out <- obj %$%
    glmgen::trendfilter(x = data.scaled$x,
                        y = data.scaled$y,
                        weights = data.scaled$weights,
                        lambda = gammas,
                        k = k,
                        thinning = thinning,
                        control = optimization.params 
    )

  obj$n.iter <- out$iter
  obj$edfs <- out$df
  obj$edf.min <- out$df[obj$i.min]
  obj$edf.1se <- out$df[obj$i.1se]
  
  gamma.pred <- case_when(
    obj$gamma.choice == "gamma.min" ~ obj$gamma.min,
    obj$gamma.choice == "gamma.1se" ~ obj$gamma.1se
  )
  
  obj$optimization.params$obj_tol <- obj$optimization.params$obj_tol * 1e-2
    
    out <- obj %$%
    glmgen::trendfilter(x = data.scaled$x,
                        y = data.scaled$y,
                        weights = data.scaled$weights,
                        lambda = gammas,
                        k = k,
                        thinning = thinning,
                        control = optimization.params 
    )
  
  obj$tf.estimate <- glmgen:::predict.trendfilter(out,
                                                  lambda = gamma.pred,
                                                  x.new = obj$x.eval / obj$x.scale
                                                  ) * obj$y.scale %>%
    as.numeric
  
  obj$fitted.values <- glmgen:::predict.trendfilter(out,
                                                    lambda = gamma.pred,
                                                    x.new = obj$data.scaled$x
                                                    ) * obj$y.scale %>%
    as.numeric

  obj$residuals <- (obj$y - obj$fitted.values)

  obj <- obj[c("x.eval","tf.estimate","validation.method","V",
               "validation.error.type","gammas","gamma.min","gamma.1se",
               "gamma.choice","errors","se.errors","edfs","edf.min","edf.1se",
               "i.min","i.1se","x","y","weights","fitted.values", "residuals",
               "k","thinning","optimization.params","n.iter","x.scale","y.scale",
               "data.scaled")]
  
  return(obj)
}


trendfilter.validate <- function(validation.index,
                                 data.folded,
                                 obj
                                 )
  {
  
  data.train <- data.folded[-validation.index] %>% bind_rows
  data.validate <- data.folded[[validation.index]]

  out <- glmgen::trendfilter(x = data.train$x,
                             y = data.train$y,
                             weights = data.train$weights,
                             k = obj$k,
                             lambda = obj$gammas,
                             thinning = obj$thinning,
                             control = obj$optimization.params
                             )
  
  tf.validate.preds <- glmgen:::predict.trendfilter(out,
                                                    lambda = obj$gammas,
                                                    x.new = data.validate$x
                                                    ) %>%
    suppressWarnings
  
  if ( obj$validation.error.type == "MSE" ){
    validation.error.mat <- (tf.validate.preds - data.validate$y) ^ 2
  }
  if ( obj$validation.error.type == "MAE" ){
    validation.error.mat <- abs(tf.validate.preds - data.validate$y)
  }
  if ( obj$validation.error.type == "WMSE" ){
    validation.error.mat <- (tf.validate.preds - data.validate$y) ^ 2 * 
    (data.validate$weights) / sum((data.validate$weights))
  }
  if ( obj$validation.error.type == "WMAE" ){
    validation.error.mat <- abs(tf.validate.preds - data.validate$y) * 
      sqrt(data.validate$weights) / sum(sqrt(data.validate$weights))
  }
  
  validation.error.sum <- colMeans(validation.error.mat) %>% as.numeric
  
  return(validation.error.sum)
}
