#' Bootstrap the optimized trend filtering estimator to obtain variability bands
#'
#' @description \loadmathjax{} \code{bootstrap.trendfilter} implements any of 
#' three possible bootstrap algorithms to obtain pointwise variability bands 
#' with a specified certainty to accompany an optimized trend filtering point 
#' estimate of a signal. 
#' @param obj An object of class 'SURE.trendfilter' or 'cv.trendfilter'.
#' @param bootstrap.method A string specifying the bootstrap method to be used. 
#' One of \code{c("nonparametric","parametric","wild")}. See Details section 
#' below for suggested use. Defaults to \code{"nonparametric"}.
#' @param alpha Specifies the width of the \code{1-alpha} pointwise variability 
#' bands. Defaults to \code{alpha = 0.05}.
#' @param B The number of bootstrap samples used to estimate the pointwise
#' variability bands. Defaults to \code{B = 250}.
#' @param full.ensemble If \code{TRUE}, the full trend filtering 
#' bootstrap ensemble is returned as an \mjeqn{n \times B}{ascii} matrix, less 
#' any columns potentially pruned post-hoc (see \code{prune} below). Defaults to 
#' \code{full.ensemble = FALSE}.
#' @param prune If \code{TRUE}, then the trend filtering bootstrap 
#' ensemble is examined for rare instances in which the optimization has 
#' stopped at zero knots (most likely erroneously), and removes them from the 
#' ensemble. Defaults to \code{TRUE}. Do not change this unless you really know 
#' what you are doing!
#' @param mc.cores Multi-core computing (for speedups): The number of cores to
#' utilize. Defaults to the number of cores detected.
#' @return An object of class 'bootstrap.trendfilter'. This is a list with the 
#' following elements:
#' \item{x.eval}{The grid of inputs the trend filtering estimate and 
#' variability bands were evaluated on.}
#' \item{tf.estimate}{The trend filtering estimate of the signal, evaluated on 
#' \code{x.eval}.}
#' \item{bootstrap.lower.perc.intervals}{Vector of lower bounds for the 
#' \code{1-alpha} pointwise variability band, evaluated on \code{x.eval}.}
#' \item{bootstrap.upper.perc.intervals}{Vector of upper bounds for the 
#' \code{1-alpha} pointwise variability band, evaluated on \code{x.eval}.}
#' \item{bootstrap.method}{The string specifying the bootstrap method that was
#' used.}
#' \item{alpha}{The 'level' of the variability bands, i.e. \code{alpha}
#' produces a \code{100*(1-alpha)}\% pointwise variability band.}
#' \item{B}{The number of bootstrap samples used to estimate the pointwise
#' variability bands.}
#' \item{tf.bootstrap.ensemble}{(Optional) The full trend filtering bootstrap 
#' ensemble as an \mjeqn{n \times B}{ascii} matrix, less any columns potentially 
#' pruned post-hoc (if \code{prune = TRUE}). If \code{full.ensemble = FALSE}, 
#' then this will return \code{NULL}.}
#' \item{df.boots}{An integer vector of the estimated number of effective 
#' degrees of freedom of each trend filtering bootstrap estimate. These should
#' all be relatively close to \code{df.min} (below).}
#' \item{prune}{If \code{TRUE}, then the trend filtering bootstrap 
#' ensemble is examined for rare instances in which the optimization has 
#' stopped at zero knots (most likely erroneously), and removes them from the 
#' ensemble.}
#' \item{n.pruned}{The number of badly-converged bootstrap trend filtering 
#' estimates pruned from the ensemble.}
#' \item{x}{The vector of the observed inputs.}
#' \item{y}{The vector of the observed outputs.}
#' \item{weights}{A vector of weights for the observed outputs. These are 
#' defined as \code{weights = 1 / sigma^2}, where \code{sigma} is a vector of 
#' standard errors of the uncertainty in the output measurements.}
#' \item{residuals}{\code{residuals = y - fitted.values}.}
#' \item{k}{The degree of the trend filtering estimator.}
#' \item{lambda}{Vector of hyperparameter values tested during validation.}
#' \item{lambda.min}{Hyperparameter value that minimizes the validation error 
#' curve.}
#' \item{df}{Integer vector of effective degrees of freedom for trend filtering
#' estimators fit during validation.}
#' \item{df.min}{The effective degrees of freedom of the optimally-tuned trend 
#' filtering estimator.}
#' \item{i.min}{The index of \code{lambda} that minimizes the validation error.}
#' \item{validation.method}{Either "SURE" or "V-fold CV".}
#' \item{error}{Vector of hyperparameter validation errors, inherited from
#' \code{obj} (either class 'SURE.trendfilter' or 'cv.trendfilter')}
#' \item{thinning}{If \code{TRUE}, then the data are 
#' preprocessed so that a smaller, better conditioned data set is used for 
#' fitting.}
#' \item{max_iter}{Maximum iterations allowed for the trend filtering 
#' convex optimization 
#' (\href{http://www.stat.cmu.edu/~ryantibs/papers/fasttf.pdf}{Ramdas and
#' Tibshirani 2016}). 
#' Increase this if the trend filtering estimate does not appear to 
#' have fully converged to a reasonable estimate of the signal.}
#' \item{obj_tol}{The tolerance used in the convex optimization stopping 
#' criterion; when the relative change in the objective function is less than 
#' this value, the algorithm terminates. Decrease this if the trend 
#' filtering estimate does not appear to have fully converged to a reasonable 
#' estimate of the signal.}
#' @details This should be a very detailed description... See
#' \href{https://academic.oup.com/mnras/article/492/3/4005/5704413}{
#' Politsch et al. (2020a)} for when each method is most appropriate.
#' @export bootstrap.trendfilter
#' @details The bootstrap method should generally be chosen according to the 
#' following criteria: \itemize{
#' \item The inputs are irregularly sampled \mjeqn{\Longrightarrow}{ascii} 
#' \code{bootstrap.method = "nonparametric"}.
#' \item The inputs are regularly sampled and the noise distribution is known 
#' \mjeqn{\Longrightarrow}{ascii} \code{bootstrap.method = "parametric"}.
#' \item The inputs are regularly sampled and the noise distribution is 
#' unknown \mjeqn{\Longrightarrow}{ascii} \code{bootstrap.method = "wild"}.}
#' @author Collin A. Politsch, \email{collinpolitsch@@gmail.com}
#' @seealso \code{\link{trendfilter}}, \code{\link{SURE.trendfilter}}, 
#' \code{\link{cv.trendfilter}}, \code{\link{relax.trendfilter}}
#' @references \enumerate{
#' \item{Politsch et al. (2020a). Trend filtering – I. A modern 
#' statistical tool for time-domain astronomy and astronomical spectroscopy. 
#' \emph{Monthly Notices of the Royal Astronomical Society}, 492(3), 
#' p. 4005-4018.
#' \href{https://academic.oup.com/mnras/article/492/3/4005/5704413}{[Link]}} \cr
#' \item{Politsch et al. (2020b). Trend Filtering – II. Denoising 
#' astronomical signals with varying degrees of smoothness. \emph{Monthly 
#' Notices of the Royal Astronomical Society}, 492(3), p. 4019-4032.
#' \href{https://academic.oup.com/mnras/article/492/3/4019/5704414}{[Link]}} \cr
#' \item{Ramdas and Tibshirani (2016). Fast and Flexible ADMM Algorithms 
#' for Trend Filtering. \emph{Journal of Computational and Graphical 
#' Statistics}, 25(3), p. 839-858.
#' \href{https://amstat.tandfonline.com/doi/abs/10.1080/10618600.2015.1054033#.XfJpNpNKju0}{[Link]}} \cr
#' \item{Arnold, Sadhanala, and Tibshirani (2014). Fast algorithms for 
#' generalized lasso problems. R package \emph{glmgen}. Version 0.0.3. 
#' \href{https://github.com/glmgen/glmgen}{[Link]}} \cr
#' \item{Hastie, Tibshirani, and Friedman (2009). The Elements of Statistical 
#' Learning: Data Mining, Inference, and Prediction. 2nd edition. Springer 
#' Series in Statistics. \href{https://web.stanford.edu/~hastie/ElemStatLearn/printings/ESLII_print12_toc.pdf}{
#' [Online print #12]}} \cr
#' \item{Mammen (1993). Bootstrap and Wild Bootstrap for High Dimensional 
#' Linear Models. \emph{The Annals of Statistics}, 21(1), p. 255-285.
#' \href{https://projecteuclid.org/journals/annals-of-statistics/volume-21/issue-1/Bootstrap-and-Wild-Bootstrap-for-High-Dimensional-Linear-Models/10.1214/aos/1176349025.full}{[Link]}} \cr
#' \item{Efron and Tibshirani (1986). Bootstrap Methods for Standard Errors, 
#' Confidence Intervals, and Other Measures of Statistical Accuracy. Statistical
#' Science, 1(1), p. 54-75.
#' \href{https://projecteuclid.org/journals/statistical-science/volume-1/issue-1/Bootstrap-Methods-for-Standard-Errors-Confidence-Intervals-and-Other-Measures/10.1214/ss/1177013815.full}{[Link]}} \cr
#' \item{Efron (1979). Bootstrap Methods: Another Look at the Jackknife.
#' \emph{The Annals of Statistics}, 7(1), p. 1-26.
#' \href{https://projecteuclid.org/journals/annals-of-statistics/volume-7/issue-1/Bootstrap-Methods-Another-Look-at-the-Jackknife/10.1214/aos/1176344552.full}{[Link]}} \cr
#' }
#' @examples 
#' #############################################################################
#' ##                    Quasar Lyman-alpha forest example                    ##
#' #############################################################################
#' 
#' # Load Lyman-alpha forest spectral observations of an SDSS quasar at redshift 
#' # z ~ 2.953. SDSS spectra are equally spaced in log10 wavelength space, 
#' # aside from some instances of masked pixels.
#' 
#' data(quasar_spec)
#' data(plotting_utilities)
#' 
#' 
#' # We are interested in denoising the observed brightness of the quasar 
#' # (measured as a 'flux' quantity) over the observed wavelength range. Since 
#' # the logarithmic wavelengths are gridded, we optimize the trend filtering 
#' # hyperparameter by minimizing the SURE estimate of fixed-input squared 
#' # prediction error. For smoothness, we use quadratic trend filtering, i.e. 
#' # the default k=2. 
#' 
#' SURE.obj <- SURE.trendfilter(x = log10.wavelength, 
#'                              y = flux, 
#'                              weights = weights)
#' 
#' 
#' # Extract the SURE error curve and optimized trend filtering estimate from 
#' # the SURE.trendfilter output
#' 
#' log.lambda <- log(SURE.obj$lambda)
#' error <- SURE.obj$error
#' log.lambda.min <- log(SURE.obj$lambda.min)
#' 
#' log10.wavelength.eval <- SURE.obj$x.eval
#' tf.estimate <- SURE.obj$tf.estimate
#' 
#' 
#' # Run a parametric bootstrap on the optimized trend filtering estimator to 
#' # obtain uncertainty bands
#' 
#' boot.out <- bootstrap.trendfilter(obj = SURE.obj, 
#'                                   bootstrap.method = "parametric")
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
#' plot(x = wavelength, y = flux, type = "l", 
#'      main = "Quasar Lyman-alpha forest", 
#'      xlab = "Observed wavelength (Angstroms)", ylab = "Flux")
#' lines(wavelength.eval, tf.estimate, col = "orange", lwd = 1.75)
#' polygon(c(wavelength.eval, rev(wavelength.eval)), 
#'         c(boot.out$bootstrap.lower.perc.intervals, 
#'         rev(boot.out$bootstrap.upper.perc.intervals)),
#'         col = transparency("orange", 90), border = NA)
#' lines(wavelength.eval, boot.out$bootstrap.lower.perc.intervals, 
#'       col = "orange", lwd = 0.5)
#' lines(wavelength.eval, boot.out$bootstrap.upper.perc.intervals, 
#'       col = "orange", lwd = 0.5)
#' legend(x = "topleft", lwd = c(1,2,8), lty = 1, cex = 0.75,
#'        col = c("black","orange", transparency("orange", 90)), 
#'        legend = c("Noisy quasar spectrum",
#'                   "Trend filtering estimate",
#'                   "95% variability band"))

#' @importFrom dplyr case_when
#' @importFrom tidyr tibble
#' @importFrom magrittr %>%
#' @importFrom parallel mclapply detectCores
#' @importFrom stats quantile
#' @importFrom glmgen trendfilter.control.list
bootstrap.trendfilter <- function(obj,
                                  bootstrap.method = c("nonparametric","parametric","wild"),
                                  alpha = 0.05, 
                                  B = 250L, 
                                  full.ensemble = FALSE,
                                  prune = TRUE,
                                  mc.cores = detectCores()
                                  )
{
  
  bootstrap.method <- match.arg(bootstrap.method)
  stopifnot( class(obj) %in% c("SURE.trendfilter", "cv.trendfilter") )
  stopifnot( alpha > 0 & alpha < 1 )
  stopifnot( B > 10 & B == round(B) )
  if ( !prune ) warning("I hope you know what you are doing!")
  obj$prune <- prune

  obj$x.scale <- mean(diff(obj$x))
  obj$y.scale <- mean(abs(obj$y)) / 10
  obj$x <- obj$x / obj$x.scale
  obj$y <- obj$y / obj$y.scale
  obj$residuals <- obj$residuals / obj$y.scale
  obj$fitted.values <- obj$fitted.values / obj$y.scale
  
  if ( is.null(obj$weights) ){
    data <- tibble(x = obj$x, y = obj$y, weights = 1,
                   fitted.values = obj$fitted.values,
                   residuals = obj$residuals
                   )
    obj$weights <- rep(1, length(obj$x))
  }else{
    obj$weights <- y.scale ^ 2 * obj$weights
    data <- tibble(x = obj$x, y = obj$y, weights = obj$weights,
                   fitted.values = obj$fitted.values,
                   residuals = obj$residuals
                   )
  }

  sampler <- case_when(
    bootstrap.method == "nonparametric" ~ list(nonparametric.resampler),
    bootstrap.method == "parametric" ~ list(parametric.sampler),
    bootstrap.method == "wild" ~ list(wild.sampler)
  )[[1]]
  
  par.func <- function(b){
    boot.tf.estimate <- tf.estimator(data = sampler(data), 
                                     obj = obj,
                                     mode = "df"
    )
    return(boot.tf.estimate)
  }
  par.out <- mclapply(1:B, par.func, mc.cores = mc.cores)
  tf.boot.ensemble <- lapply(X = 1:B, 
                             FUN = function(X) par.out[[X]][["tf.estimate"]]) %>%
    unlist %>% 
    matrix(nrow = length(obj$x.eval)) 
  
  obj$df.boots <- lapply(X = 1:B, FUN = function(X) par.out[[X]][["df"]]) %>%
    unlist %>%
    as.integer
  obj$n.pruned <- (B - ncol(tf.boot.ensemble)) %>% as.integer
  obj$bootstrap.lower.perc.intervals <- apply(tf.boot.ensemble, 1, quantile, 
                                              probs = alpha / 2)
  obj$bootstrap.upper.perc.intervals <- apply(tf.boot.ensemble, 1, quantile, 
                                              probs = 1 - alpha / 2)
  obj <- c(obj, list(bootstrap.method = bootstrap.method, alpha = alpha, B = B))
  
  obj$x.eval <- obj$x.eval * obj$x.scale
  obj$tf.estimate <- obj$tf.estimate * obj$y.scale
  obj$x <- obj$x * obj$x.scale
  obj$y <- obj$y * obj$y.scale
  obj$weights <- obj$weights / obj$y.scale ^ 2
  obj$fitted.values <- obj$fitted.values * obj$y.scale
  obj$residuals <- (obj$y - obj$fitted.values) * obj$y.scale
  
  if ( full.ensemble ){
    obj$tf.bootstrap.ensemble <- tf.boot.ensemble
  }else{
    obj <- c(obj, list(tf.bootstrap.ensemble = NULL))
  }
  
  obj <- obj[c("x.eval","tf.estimate","bootstrap.lower.perc.intervals",
               "bootstrap.upper.perc.intervals","bootstrap.method","alpha","B",
               "df.boots","tf.bootstrap.ensemble","prune","n.pruned","x","y",
               "weights","fitted.values","residuals","k","lambda","lambda.min",
               "df","df.min","i.min","validation.method","error","thinning",
               "max_iter","obj_tol")
             ]
  class(obj) <- "bootstrap.trendfilter"
  
  return(obj)
}


#' @importFrom dplyr %>%
tf.estimator <- function(data, 
                         obj,
                         mode = "lambda"
                         )
  {
  
  optimization.controls <- glmgen::trendfilter.control.list(max_iter = obj$max_iter,
                                                            obj_tol = obj$obj_tol
                                                            )

  if ( mode == "df" ){
    tf.fit <- glmgen::trendfilter(x = data$x,
                                  y = data$y,
                                  weights = data$weights,
                                  k = obj$k,
                                  lambda = obj$lambda,
                                  thinning = obj$thinning,
                                  control = optimization.controls
                                  )
    
    i.min <- which.min( abs(tf.fit$df - obj$df.min) )
    lambda.min <- obj$lambda[i.min]
    df.min <- tf.fit$df[i.min]
    
    if ( obj$prune & df.min <= 2 ){
      return(list(tf.estimate = integer(0), df = NA))
    }
    
  }else{
    tf.fit <- glmgen::trendfilter(x = data$x,
                                  y = data$y,
                                  weights = data$weights,
                                  k = obj$k,
                                  lambda = obj$lambda.min,
                                  thinning = obj$thinning,
                                  control = optimization.controls
                                  )
    
    lambda.min <- obj$lambda.min
    df.min <- tf.fit$df %>% as.integer
  }

  tf.estimate <- glmgen:::predict.trendfilter(object = tf.fit, 
                                              x.new = obj$x.eval / obj$x.scale, 
                                              lambda = lambda.min
                                              ) %>% as.numeric
  
  return(list(tf.estimate = tf.estimate * obj$y.scale, df = df.min))
}

#' @importFrom dplyr slice_sample
nonparametric.resampler <- function(data){
  slice_sample(data, n = nrow(data), replace = TRUE)
}


#' @importFrom stats rnorm
#' @importFrom dplyr %>% mutate n
parametric.sampler <- function(data){
  data %>% mutate(y = fitted.values + rnorm(n = n(), sd = 1 / sqrt(weights)))
}


#' @importFrom dplyr %>% mutate n
wild.sampler <- function(data){
  data %>% mutate(y = fitted.values + residuals *
                    sample(x = c(
                      (1 + sqrt(5)) / 2, 
                      (1 - sqrt(5)) / 2
                      ), 
                           size = n(), replace = TRUE,
                           prob = c(
                             (1 + sqrt(5)) / (2 * sqrt(5)),
                             (sqrt(5) - 1) / (2 * sqrt(5))
                                    )
                           )
                  )
}
