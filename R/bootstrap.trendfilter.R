#' Bootstrap the optimized trend filtering estimator to obtain variability bands
#'
#' @description \loadmathjax{} \code{bootstrap.trendfilter} implements any of 
#' three possible bootstrap algorithms to obtain pointwise variability bands 
#' with a specified certainty to accompany an optimized trend filtering point 
#' estimate of a signal. 
#' @param obj An object of class '\link{SURE.trendfilter}' or 
#' '\link{cv.trendfilter}'.
#' @param bootstrap.algorithm A string specifying the bootstrap algorithm to be 
#' used. One of \code{c("nonparametric","parametric","wild")}. See Details 
#' section below for suggested use. Defaults to \code{"nonparametric"}.
#' @param alpha Specifies the width of the \code{1-alpha} pointwise variability 
#' bands. Defaults to \code{alpha = 0.05}.
#' @param B The number of bootstrap samples used to estimate the pointwise
#' variability bands. Defaults to \code{B = 100}. Increase this for more precise
#' bands (e.g. the final analysis you intend to publish).
#' @param full.ensemble logical. If \code{TRUE}, the full trend filtering 
#' bootstrap ensemble is returned as an \mjeqn{n \times B}{ascii} matrix, less 
#' any columns potentially pruned post-hoc (see \code{prune} below). Defaults to 
#' \code{full.ensemble = FALSE}.
#' @param prune logical. If \code{TRUE}, then the trend filtering bootstrap 
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
#' \item{bootstrap.lower.band}{Vector of lower bounds for the 
#' \code{1-alpha} pointwise variability band, evaluated on \code{x.eval}.}
#' \item{bootstrap.upper.band}{Vector of upper bounds for the 
#' \code{1-alpha} pointwise variability band, evaluated on \code{x.eval}.}
#' \item{bootstrap.algorithm}{The string specifying the bootstrap algorithms 
#' that was used.}
#' \item{alpha}{The 'level' of the variability bands, i.e. \code{alpha}
#' produces a \code{100*(1-alpha)}\% pointwise variability band.}
#' \item{B}{The number of bootstrap samples used to estimate the pointwise
#' variability bands.}
#' \item{tf.bootstrap.ensemble}{(Optional) If \code{full.ensemble = TRUE}, the 
#' full trend filtering bootstrap ensemble as an \mjeqn{n \times B}{ascii} 
#' matrix, less any columns potentially pruned post-hoc 
#' (if \code{prune = TRUE}). If \code{full.ensemble = FALSE}, then this will 
#' return \code{NULL}.}
#' \item{edf.boots}{An integer vector of the estimated number of effective 
#' degrees of freedom of each trend filtering bootstrap estimate. These should
#' all be relatively close to \code{edf.min} (below).}
#' \item{prune}{logical. If \code{TRUE}, then the trend filtering bootstrap 
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
#' \item{gammas}{Vector of hyperparameter values tested during validation.}
#' \item{gammas.min}{Hyperparameter value that minimizes the validation error 
#' curve.}
#' \item{edf}{Integer vector of effective degrees of freedom for trend filtering
#' estimators fit during validation.}
#' \item{edf.min}{The effective degrees of freedom of the optimally-tuned trend 
#' filtering estimator.}
#' \item{i.min}{The index of \code{gammas} that minimizes the validation error.}
#' \item{validation.method}{Either "SURE" or "V-fold CV".}
#' \item{error}{Vector of hyperparameter validation errors, inherited from
#' \code{obj} (either class 'SURE.trendfilter' or 'cv.trendfilter')}
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
#' @details This will contain a very detailed description... See
#' \href{https://academic.oup.com/mnras/article/492/3/4005/5704413}{
#' Politsch et al. (2020a)} for the full algorithms. The bootstrap algorithm 
#' should generally be chosen according to the following criteria: \itemize{
#' \item The inputs are irregularly sampled \mjeqn{\Longrightarrow}{ascii} 
#' \code{bootstrap.algorithm = "nonparametric"}.
#' \item The inputs are regularly sampled and the noise distribution is known 
#' \mjeqn{\Longrightarrow}{ascii} \code{bootstrap.algorithm = "parametric"}.
#' \item The inputs are regularly sampled and the noise distribution is 
#' unknown \mjeqn{\Longrightarrow}{ascii} \code{bootstrap.algorithm = "wild"}.}
#' @export bootstrap.trendfilter
#' @author Collin A. Politsch, \email{collinpolitsch@@gmail.com}
#' @seealso \code{\link{SURE.trendfilter}}, \code{\link{cv.trendfilter}}
#' @references 
#' \strong{Trend filtering and the various bootstraps in practice}
#' \enumerate{
#' \item{Politsch et al. (2020a). Trend filtering – I. A modern 
#' statistical tool for time-domain astronomy and astronomical spectroscopy. 
#' \emph{Monthly Notices of the Royal Astronomical Society}, 492(3), 
#' p. 4005-4018.
#' \href{https://academic.oup.com/mnras/article/492/3/4005/5704413}{[Link]}} \cr
#' \item{Politsch et al. (2020b). Trend Filtering – II. Denoising 
#' astronomical signals with varying degrees of smoothness. \emph{Monthly 
#' Notices of the Royal Astronomical Society}, 492(3), p. 4019-4032.
#' \href{https://academic.oup.com/mnras/article/492/3/4019/5704414}{[Link]}} \cr
#' }
#' \strong{The Bootstrap (and variations)}
#' \enumerate{
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
#' \strong{Trend filtering optimization algorithm}
#' \enumerate{
#' \item{Ramdas and Tibshirani (2016). Fast and Flexible ADMM Algorithms 
#' for Trend Filtering. \emph{Journal of Computational and Graphical 
#' Statistics}, 25(3), p. 839-858.
#' \href{https://amstat.tandfonline.com/doi/abs/10.1080/10618600.2015.1054033#.XfJpNpNKju0}{[Link]}} \cr
#' \item{Arnold, Sadhanala, and Tibshirani (2014). Fast algorithms for 
#' generalized lasso problems. R package \emph{glmgen}. Version 0.0.3. 
#' \href{https://github.com/glmgen/glmgen}{[Link]}}
#' (Implementation of Ramdas and Tibshirani algorithm) \cr
#' }
#' @examples 
#' #############################################################################
#' ##                    Quasar Lyman-alpha forest example                    ##
#' #############################################################################
#' 
#' # A quasar is an extremely luminous galaxy with an active supermassive black 
#' # hole at its center. Absorptions in the spectra of quasars at vast 
#' # cosmological distances from our galaxy reveal the presence of a gaseous 
#' # medium permeating the entirety of intergalactic space -- appropriately 
#' # named the 'intergalactic medium'. These absorptions allow astronomers to 
#' # study the structure of the Universe using the distribution of these 
#' # absorptions in quasar spectra. Particularly important is the 'forest' of 
#' # absorptions that arise from the Lyman-alpha spectral line, which traces 
#' # the presence of electrically neutral hydrogen in the intergalactic medium.
#' #
#' # Here, we are interested in denoising the Lyman-alpha forest of a quasar 
#' # spectroscopically measured by the Sloan Digital Sky Survey. SDSS spectra 
#' # are equally spaced in log10 wavelength space, aside from some instances of 
#' # masked pixels.
#' 
#' data(quasar_spec)
#'
#' # head(data)
#' #
#' # | log10.wavelength|       flux|   weights|
#' # |----------------:|----------:|---------:|
#' # |           3.5529|  0.4235348| 0.0417015|
#' # |           3.5530| -2.1143005| 0.1247811|
#' # |           3.5531| -3.7832341| 0.1284383|
#' 
#' SURE.obj <- SURE.trendfilter(x = data$log10.wavelength, 
#'                              y = data$flux, 
#'                              weights = data$weights)
#' 
#' 
#' # Extract the estimated hyperparameter error curve and optimized trend 
#' # filtering estimate from the `SURE.trendfilter` output, and transform the 
#' # input grid to wavelength space (in Angstroms).
#' 
#' log.gammas <- log(SURE.obj$gammas)
#' errors <- SURE.obj$errors
#' log.gamma.min <- log(SURE.obj$gamma.min)
#' 
#' wavelength <- 10 ^ (SURE.obj$x)
#' wavelength.eval <- 10 ^ (SURE.obj$x.eval)
#' tf.estimate <- SURE.obj$tf.estimate
#'
#' 
#' # Run a parametric bootstrap on the optimized trend filtering estimator to 
#' # obtain uncertainty bands
#' 
#' boot.out <- bootstrap.trendfilter(obj = SURE.obj, 
#'                                   bootstrap.algorithm = "parametric")
#' 
#' 
#' # Plot the results
#'
#' par(mfrow = c(2,1), mar = c(5,4,2.5,1) + 0.1)
#' plot(x = log.gammas, y = errors, main = "SURE error curve", 
#'      xlab = "log(gamma)", ylab = "SURE error")
#' abline(v = log.gamma.min, lty = 2, col = "blue3")
#' text(x = log.gamma.min, y = par("usr")[4], 
#'      labels = "optimal gamma", pos = 1, col = "blue3")
#' 
#' plot(x = wavelength, y = SURE.obj$y, type = "l", 
#'      main = "Quasar Lyman-alpha forest", 
#'      xlab = "Observed wavelength (Angstroms)", ylab = "Flux")
#' polygon(c(wavelength.eval, rev(wavelength.eval)), 
#'         c(boot.out$bootstrap.lower.band, 
#'         rev(boot.out$bootstrap.upper.band)),
#'         col = transparency("orange", 90), border = NA)
#' lines(wavelength.eval, boot.out$bootstrap.lower.band, 
#'       col = "orange", lwd = 0.5)
#' lines(wavelength.eval, boot.out$bootstrap.upper.band, 
#'       col = "orange", lwd = 0.5)
#' lines(wavelength.eval, tf.estimate, col = "orange", lwd = 2.5)
#' legend(x = "topleft", lwd = c(1,2,8), lty = 1, cex = 0.75,
#'        col = c("black","orange", transparency("orange", 90)), 
#'        legend = c("Noisy quasar spectrum",
#'                   "Trend filtering estimate",
#'                   "95% variability band"))

#' @importFrom dplyr case_when mutate slice_sample n
#' @importFrom tidyr tibble
#' @importFrom magrittr %>%
#' @importFrom parallel mclapply detectCores
#' @importFrom stats quantile rnorm
bootstrap.trendfilter <- function(obj,
                                  bootstrap.algorithm = c("nonparametric","parametric","wild"),
                                  alpha = 0.05, 
                                  B = 100L, 
                                  full.ensemble = FALSE,
                                  prune = TRUE,
                                  mc.cores = detectCores()
                                  )
{
  
  stopifnot( class(obj) %in% c("SURE.trendfilter", "cv.trendfilter") )
  bootstrap.algorithm <- match.arg(bootstrap.algorithm)
  stopifnot( alpha > 0 & alpha < 1 )
  stopifnot( B > 20 & B == round(B) )
  if ( !prune ) warning("I hope you know what you are doing!")
  if ( mc.cores != round(mc.cores) ) stop("mc.cores must be a positive integer.")
  if ( mc.cores > detectCores() ){
    warning(paste0("Your machine only has ", detectCores(), " cores. Adjusting mc.cores accordingly."))
    mc.cores <- detectCores()
  }

  sampler <- case_when(
    bootstrap.algorithm == "nonparametric" ~ list(nonparametric.resampler),
    bootstrap.algorithm == "parametric" ~ list(parametric.sampler),
    bootstrap.algorithm == "wild" ~ list(wild.sampler)
  )[[1]]
  
  obj$prune <- prune
  data.scaled <- obj$data.scaled
  
  par.func <- function(b){
    boot.tf.estimate <- tf.estimator(data = sampler(data.scaled), 
                                     obj = obj,
                                     mode = "edf"
    )
    return(boot.tf.estimate)
  }
  par.out <- mclapply(1:B, par.func, mc.cores = mc.cores)
  tf.boot.ensemble <- lapply(X = 1:B, 
                             FUN = function(X) par.out[[X]][["tf.estimate"]]) %>%
    unlist %>% 
    matrix(nrow = length(obj$x.eval)) 
  
  obj$edf.boots <- lapply(X = 1:B, FUN = function(X) par.out[[X]][["edf"]]) %>%
    unlist %>%
    as.integer
  obj$n.pruned <- (B - ncol(tf.boot.ensemble)) %>% as.integer
  obj$bootstrap.lower.band <- apply(tf.boot.ensemble, 1, quantile, 
                                              probs = alpha / 2) 
  obj$bootstrap.upper.band <- apply(tf.boot.ensemble, 1, quantile, 
                                              probs = 1 - alpha / 2)
  obj <- c(obj, list(bootstrap.algorithm = bootstrap.algorithm, alpha = alpha, B = B))
  
  
  if ( full.ensemble ){
    obj$tf.bootstrap.ensemble <- tf.boot.ensemble
  }else{
    obj <- c(obj, list(tf.bootstrap.ensemble = NULL))
  }
  
  obj <- obj[c("x.eval","tf.estimate","bootstrap.lower.band",
               "bootstrap.upper.band","bootstrap.algorithm","alpha","B",
               "edf.boots","tf.bootstrap.ensemble","prune","n.pruned","x","y",
               "weights","fitted.values","residuals","k","gammas","gamma.min",
               "edfs","edf.min","i.min","validation.method","errors","thinning",
               "optimization.params","n.iter","x.scale","y.scale","data.scaled")
             ]
  class(obj) <- "bootstrap.trendfilter"
  
  return(obj)
}


tf.estimator <- function(data, 
                         obj,
                         mode = "gamma"
                         )
  {
  
  optimization.controls <- glmgen::trendfilter.control.list(max_iter = obj$max_iter,
                                                            obj_tol = obj$obj_tol
                                                            )

  if ( mode == "edf" ){
    tf.fit <- glmgen::trendfilter(x = data$x,
                                  y = data$y,
                                  weights = data$weights,
                                  k = obj$k,
                                  lambda = obj$gammas,
                                  thinning = obj$thinning,
                                  control = optimization.controls
                                  )
    
    i.min <- which.min( abs(tf.fit$df - obj$edf.min) )
    gammas.min <- obj$gammas[i.min]
    edf.min <- tf.fit$df[i.min]
    
    if ( obj$prune & df.min <= 2 ){
      return(list(tf.estimate = integer(0), df = NA))
    }
    
  }else{
    tf.fit <- glmgen::trendfilter(x = data$x,
                                  y = data$y,
                                  weights = data$weights,
                                  k = obj$k,
                                  lambda = obj$gamma.min,
                                  thinning = obj$thinning,
                                  control = optimization.controls
                                  )
    
    gamma.min <- obj$gamma.min
    edf.min <- tf.fit$df
  }

  tf.estimate <- glmgen:::predict.trendfilter(object = tf.fit, 
                                              x.new = obj$x.eval / obj$x.scale, 
                                              lambda = gamma.min
                                              ) %>% as.numeric
  
  return(list(tf.estimate = tf.estimate * obj$y.scale, edf = edf.min))
}


nonparametric.resampler <- function(data){
  slice_sample(data, n = nrow(data), replace = TRUE)
}


parametric.sampler <- function(data){
  data %>% mutate(y = fitted.values + rnorm(n = n(), sd = 1 / sqrt(weights)))
}


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
