#' Bootstrap the optimized trend filtering estimator to obtain variability bands
#'
#' @description \code{bootstrap.trendfilter} implements any of three possible 
#' bootstrap algorithms to obtain pointwise variability bands with a specified 
#' certainty to accompany an optimized trend filtering point estimate of a 
#' signal. 
#' @param obj An object of class 'SURE.trendfilter' or 'cv.trendfilter'.
#' @param bootstrap.method A string specifying the bootstrap method to be used. 
#' One of \code{c("nonparametric","parametric","wild")}. See Details section 
#' below for suggested use. Defaults to \code{"nonparametric"}.
#' @param alpha Specifies the width of the \code{1-alpha} pointwise variability 
#' bands. Defaults to \code{alpha = 0.05}.
#' @param B The number of bootstrap samples used to estimate the pointwise
#' variability bands. Defaults to \code{B = 250}.
#' @param full.ensemble (logical) If \code{TRUE}, the full trend filtering 
#' bootstrap ensemble is returned as an \eqn{n x B} matrix. Defaults to 
#' \code{FALSE}.
#' @param prune (logical) If \code{TRUE}, then the trend filtering bootstrap 
#' ensemble is examined for rare instances in which the optimization has stopped 
#' at zero knots (most likely erroneously), and removes them from the ensemble. 
#' Defaults to \code{TRUE}. Do not change this unless you really know what you 
#' are doing.
#' @param mc.cores Multi-core computing (for speedups): The number of cores to
#' utilize. If 4 or more cores are detected, then the default is to utilize
#' \code{available.cores - 2}. Else, \code{mc.cores = 1}.
#' @return An object of class 'bootstrap.trendfilter'. This is a list with the 
#' following elements:
#' \item{x.eval}{The grid of inputs the trend filtering estimate and variability 
#' bands were evaluated on.}
#' \item{tf.estimate}{The trend filtering estimate of the signal, evaluated on 
#' \code{x.eval}.}
#' \item{bootstrap.lower.perc.intervals}{Vector of lower bounds for the 1-alpha 
#' pointwise variability band, evaluated on \code{x.eval}.}
#' \item{bootstrap.upper.perc.intervals}{Vector of upper bounds for the 1-alpha 
#' pointwise variability band, evaluated on \code{x.eval}.}
#' \item{bootstrap.method}{The string specifying the bootstrap method that was
#' used.}
#' \item{alpha}{The 'level' of the variability bands, i.e. \code{alpha}
#' produces a \eqn{100*(1-\alpha)}\% pointwise variability band.}
#' \item{B}{The number of bootstrap samples used to estimate the pointwise
#' variability bands.}
#' \item{tf.bootstrap.ensemble}{(Optional) The full trend filtering bootstrap 
#' ensemble as an \eqn{n x B} matrix. If \code{full.ensemble = FALSE}, then 
#' this will return \code{NULL}.}
#' \item{prune}{(logical) If \code{TRUE}, then the trend filtering bootstrap 
#' ensemble is examined for rare instances in which the optimization has 
#' stopped at zero knots (most likely erroneously), and removes them from the 
#' ensemble.}
#' \item{n.pruned}{The number of badly-converged bootstrap trend filtering 
#' estimates pruned from the ensemble.}
#' \item{x}{The vector of the observed inputs.}
#' \item{y}{The vector of the observed outputs.}
#' \item{weights}{A vector of weights for the observed outputs. These are
#' defined as \code{weights[i] = 1 / sigma[i]^2}, where \code{sigma} is a vector 
#' of standard errors of the uncertainty in the measured outputs.}
#' \item{fitted.values}{The trend filtering estimate of the signal, evaluated at
#' the observed inputs \code{x}.}
#' \item{residuals}{\code{residuals = y - fitted.values}.}
#' \item{k}{(Integer) The degree of the trend filtering estimator.}
#' \item{lambda}{Vector of hyperparameter values tested during validation.}
#' \item{lambda.min}{Hyperparameter value that minimizes the validation error 
#' curve.}
#' \item{df}{Vector of effective degrees of freedom for trend filtering
#' estimators fit during validation.}
#' \item{df.min}{The effective degrees of freedom of the optimally-tuned trend 
#' filtering estimator.}
#' \item{i.min}{The index of \code{lambda} that minimizes the validation error.}
#' \item{validation.method}{Either "SURE" or "cv".}
#' \item{error}{Vector of hyperparameter validation errors, inherited from
#' \code{obj} (either class 'SURE.trendfilter' or 'cv.trendfilter')}
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
#' @export bootstrap.trendfilter
#' @details The bootstrap method should generally be chosen according to the 
#' following criteria: \itemize{
#' \item The inputs are irregularly sampled -> 
#' \code{bootstrap.method = "nonparametric"}.
#' \item The inputs are regularly sampled and the noise distribution is known –> 
#' \code{bootstrap.method = "parametric"}.
#' \item The inputs are regularly sampled and the noise distribution is 
#' unknown –> \code{bootstrap.method = "wild"}.}
#' See \href{https://academic.oup.com/mnras/article/492/3/4005/5704413}{
#' Politsch et al. (2020)} for more details.
#' @author Collin A. Politsch, \email{collinpolitsch@@gmail.com}
#' @seealso \code{\link{SURE.trendfilter}}, \code{\link{cv.trendfilter}}
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
#' # Run a parametric bootstrap the optimized trend filtering estimator to 
#' # obtain uncertainty bands
#' 
#' boot.out <- bootstrap.trendfilter(obj = SURE.obj, bootstrap.method = "parametric")
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
#' polygon(c(wavelength.eval, rev(wavelength.eval)), 
#'         c(boot.out$bootstrap.lower.perc.intervals, 
#'         rev(boot.out$bootstrap.upper.perc.intervals)),
#'         col = transparency("orange", 90), border=NA)
#' lines(wavelength.eval, boot.out$bootstrap.lower.perc.intervals, 
#'       col = "orange", lwd = 0.5)
#' lines(wavelength.eval, boot.out$bootstrap.upper.perc.intervals, 
#'       col = "orange", lwd = 0.5)
#' legend(x = "topleft", lwd = c(1,2,8), lty = 1, cex = 0.75,
#'        col = c("black","orange", transparency("orange", 90)), 
#'        legend = c("Noisy quasar spectrum",
#'                   "Trend filtering estimate",
#'                   "95 percent variability band"))

#' @importFrom stats quantile
bootstrap.trendfilter <- function(obj,
                                  bootstrap.method = c("nonparametric","parametric","wild"),
                                  alpha = 0.05, 
                                  B = 250L, 
                                  full.ensemble = FALSE,
                                  prune = TRUE,
                                  mc.cores = max(c(parallel::detectCores() - 2), 1)
                                  )
{
  
  bootstrap.method <- match.arg(bootstrap.method)
  stopifnot( class(obj) %in% c("SURE.trendfilter", "cv.trendfilter") )
  stopifnot( alpha > 0 & alpha < 1 )
  stopifnot( B > 10 & B == round(B) )
  if ( !prune ) warning("I hope you know what you are doing!")
  obj$prune <- prune

  if ( is.null(obj$weights) ){
    data <- data.frame(x = obj$x, y = obj$y, weights = 1,
                       fitted.values = obj$fitted.values, residuals = obj$residuals
                       )
    obj$weights <- rep(1, length(obj$x))
  }else{
    data <- data.frame(x = obj$x, y = obj$y, weights = obj$weights,
                       fitted.values = obj$fitted.values, residuals = obj$residuals
                       )
  }

  sampler <- case_when(
    bootstrap.method == "nonparametric" ~ list(nonparametric.resampler),
    bootstrap.method == "parametric" ~ list(parametric.sampler),
    bootstrap.method == "wild" ~ list(wild.sampler)
  )[[1]]
  
  if ( mc.cores == 1 ){
    tf.boot.ensemble <- matrix(unlist(replicate(B,
                                                tf.estimator(data = sampler(data), 
                                                             obj = obj,
                                                             mode = "df"
                                                             )
                                                )
                                      ),
                             nrow = length(obj$x.eval)
                             )
  }else{
    par.func <- function(b){
      boot.tf.estimate <- tf.estimator(data = sampler(data), 
                                       obj = obj,
                                       mode = "df"
      )
      return(boot.tf.estimate)
    }
    tf.boot.ensemble <- matrix(unlist(parallel::mclapply(1:B, par.func, mc.cores = mc.cores)), 
                               nrow = length(obj$x.eval)
                               )
  }
  
  obj$n.pruned <- B - ncol(tf.boot.ensemble)
  obj$bootstrap.lower.perc.intervals <- apply(tf.boot.ensemble, 1, quantile, probs = alpha/2)
  obj$bootstrap.upper.perc.intervals <- apply(tf.boot.ensemble, 1, quantile, probs = 1-alpha/2)
  obj <- c(obj, list(bootstrap.method = bootstrap.method, alpha = alpha, B = B) )
  
  if ( full.ensemble ){
    obj$tf.bootstrap.ensemble <- tf.boot.ensemble
  }else{
    obj <- c(obj, list(tf.bootstrap.ensemble = NULL))
  }
  
  obj <- obj[c("x.eval","tf.estimate","bootstrap.lower.perc.intervals","bootstrap.upper.perc.intervals",
               "bootstrap.method","alpha","B","tf.bootstrap.ensemble","prune","n.pruned","x","y",
               "weights", "fitted.values","residuals","k","lambda","lambda.min","df","df.min","i.min",
               "validation.method","error","thinning","max_iter","obj_tol")
             ]
  class(obj) <- "bootstrap.trendfilter"
  
  return(obj)
}


tf.estimator <- function(data, 
                         obj,
                         mode = "lambda"
                         )
  {

  if ( mode == "df" ){
    tf.fit <- trendfilter(x = data$x, 
                          y = data$y, 
                          weights = data$weights, 
                          k = obj$k,
                          lambda = obj$lambda, 
                          thinning = obj$thinning,
                          control = trendfilter.control.list(max_iter = obj$max_iter,
                                                             obj_tol = obj$obj_tol
                                                             )
                          )
    
    i.min <- which.min( abs(tf.fit$df - obj$df.min) )
    lambda.min <- obj$lambda[i.min]
    
    if ( obj$prune & obj$df[i.min] <= 2 ){
      return(NULL)
    }
    
  }else{
    tf.fit <- trendfilter(x = data$x, 
                          y = data$y, 
                          weights = data$weights,
                          k = obj$k,
                          lambda = obj$lambda.min, 
                          thinning = obj$thinning,
                          control = trendfilter.control.list(max_iter = obj$max_iter,
                                                             obj_tol = obj$obj_tol
                                                             )
                          )
    
    lambda.min <- obj$lambda.min
  }

  tf.estimate <- glmgen:::predict.trendfilter(object = tf.fit, 
                                              x.new = obj$x.eval, 
                                              lambda = lambda.min
                                              ) %>% as.numeric
  
  return(tf.estimate)
}


nonparametric.resampler <- function(data){
  resampled.data <- slice_sample(data, n = nrow(data), replace = TRUE)
  return(resampled.data)
}


#' @importFrom stats rnorm
parametric.sampler <- function(data){
  boot.sample <- data$fitted.values + rnorm(nrow(data), sd = 1 / sqrt(data$weights))
  return(data.frame(x = data$x, y = boot.sample, weights = data$weights))
}


wild.sampler <- function(data){
  wild.boot.residuals <- data$residuals * sample(x = c((1+sqrt(5))/2, 1-sqrt(5)/2), 
                                                 size = nrow(data), 
                                                 replace = TRUE, 
                                                 prob = c((1+sqrt(5))/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5)))
                                                 )
  wild.boot.sample <- data$fitted.values + wild.boot.residuals
  return(data.frame(x = data$x, y = wild.boot.sample, weights = data$weights))
}
