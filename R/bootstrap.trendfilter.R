#' Bootstrap the optimized trend filtering estimator to obtain variability bands
#'
#' @description \code{bootstrap.trendfilter} implements any of three possible 
#' bootstrap algorithms to obtain pointwise variability bands with a specified 
#' certainty to accompany an optimized trend filtering point estimate of a 
#' signal. 
#' @param obj An object of class \code{SURE.trendfilter} or 
#' \code{cv.trendfilter}.
#' @param x.eval.grid Grid of inputs to evaluate the variability bands on. 
#' Defaults to the observed inputs.
#' @param bootstrap.method A string specifying the bootstrap method to be used. 
#' See Details section below for suggested use. Defaults to 
#' \code{bootstrap.method = "nonparametric"}.
#' @param alpha Specifies the width of the \code{1-alpha} pointwise variability 
#' bands. Defaults to \code{alpha = 0.05}.
#' @param B The number of bootstrap samples used to estimate the pointwise
#' variability bands. Defaults to \code{B = 1000}.
#' @param full.ensemble Return the full bootstrap ensemble as an \code{n x B} 
#' matrix. Defaults to \code{full.ensemble = FALSE}.
#' @param mc.cores Multi-core computing (for speedups): The number of cores to
#' utilize. If 4 or more cores are detected, then the default is to utilize
#' \code{n.cores - 2}. Else, \code{mc.cores = 1}.
#' @return A list with the following elements:
#' \item{bootstrap.lower.perc.intervals}{Vector of lower bounds for the 1-alpha 
#' pointwise variability band.}
#' \item{bootstrap.upper.perc.intervals}{Vector of upper bounds for the 1-alpha 
#' pointwise variability band.}
#' \item{tf.boot.ensemble}{(Optional) The full bootstrap ensemble as an 
#' \code{n x B} matrix.}
#' @export bootstrap.trendfilter
#' @details The bootstrap method should generally be chosen according to the 
#' following criteria: \itemize{
#' \item The inputs are irregularly sampled –> 
#' \code{bootstrap.method = "nonparametric"}
#' \item The inputs are regularly sampled and the noise distribution is known –> 
#' \code{bootstrap.method = "parametric"} 
#' \item The inputs are regularly sampled and the noise distribution is 
#' unknown –> \code{bootstrap.method = "wild"}}
#' @author Collin A. Politsch, \email{collinpolitsch@@gmail.com}
#' @seealso \code{\link{SURE.trendfilter}}
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
#' lambda.grid <- exp(seq(-10, 5, length = 200))
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
#' 
#' boot.out <- bootstrap.trendfilter(obj = SURE.obj,
#'                                   bootstrap.method = "parametric",
#'                                   x.eval.grid = x.eval.grid
#'                                   )
#'                                   
#'                                   
#' # Superpose a transparent 95% variability 'envelope' on top of the trend 
#' # filtering point estimate
#' 
#' # transparency() adds transparency to a color. Define transparency with an 
#' # integer between 0 and 255, 0 being fully transparent and 255 being fully 
#' # visible.
#' 
#' transparency <- function(color, trans){
#' 
#'     num2hex <- function(x)
#'     {
#'       hex <- unlist(strsplit("0123456789ABCDEF",split=""))
#'       return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
#'     }
#'     rgb <- rbind(col2rgb(color),trans)
#'     res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
#'     return(res)
#' }
#'                                   
#' polygon(c(wavelength.eval.grid, rev(wavelength.eval.grid)), 
#'         c(boot.out$bootstrap.lower.perc.intervals, 
#'         rev(boot.out$bootstrap.upper.perc.intervals)),
#'         col = transparency("orange", 90), border=NA)
#' lines(wavelength.eval.grid, boot.out$bootstrap.lower.perc.intervals, 
#'       col = "orange", lwd = 0.5)
#' lines(wavelength.eval.grid, boot.out$bootstrap.upper.perc.intervals, 
#'       col = "orange", lwd = 0.5)
#' legend(x = "topleft", lwd = c(2,8), lty = 1, 
#'        col = c("orange", transparency("orange", 90)), 
#'        legend = c("Trend filtering estimate", "95 percent variability band"))

#' @importFrom stats quantile
#' @importFrom dplyr case_when
bootstrap.trendfilter <- function(obj,
                                  bootstrap.method = "nonparametric", 
                                  x.eval.grid = NULL, 
                                  alpha = 0.05, 
                                  B = 1000L, 
                                  full.ensemble = FALSE,
                                  mc.cores = max(c(parallel::detectCores() - 2), 1)
                                  )
{
  
  if ( !(class(obj) %in% c("SURE.trendfilter", "cv.trendfilter")) ){
    stop("obj must be an object of class 'SURE.trendfilter' or 'cv.trendfilter'.")
  }
  
  if ( is.null(obj$weights) ){
    data <- data.frame(x = obj$x, y = obj$y, weights = 1)
    obj$weights <- rep(1, length(obj$x))
  }else{
    data <- data.frame(x = obj$x, y = obj$y, weights = obj$weights)
  }

  if ( is.null(x.eval.grid) ){
    x.eval.grid <- seq(min(obj$x), max(obj$x), length = 1500)
  }
  
  obj$x.eval.grid <- x.eval.grid
  data$tf.estimate <- tf.estimator(data = data, 
                                   obj = obj,
                                   mode = "lambda",
                                   x.eval.grid = data$x
                                   )
  data$residuals <- data$y - data$tf.estimate
  obj$tf.estimate <- data$tf.estimate
  obj$residuals <- data$residuals
  
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
                             nrow = length(obj$x.eval.grid)
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
                               nrow = length(obj$x.eval.grid)
                               )
  }
  
  obj$bootstrap.lower.perc.intervals <- apply(tf.boot.ensemble, 1, quantile, probs = alpha/2)
  obj$bootstrap.upper.perc.intervals <- apply(tf.boot.ensemble, 1, quantile, probs = 1-alpha/2)
  obj <- c(obj, list(bootstrap.method = bootstrap.method,
                     alpha = alpha, 
                     B = B
                     )
           )
  
  if ( full.ensemble ){
    obj$tf.bootstrap.ensemble <- tf.boot.ensemble
  }else{
    obj$tf.bootstrap.ensemble <- NULL
  }
  
  class(obj) <- "bootstrap.trendfilter"
  
  return(obj)
}

#' @importFrom stats predict
tf.estimator <- function(data, 
                         obj = obj,
                         mode = "lambda",
                         x.eval.grid = NULL)
  {
  
  if ( is.null(x.eval.grid) ){
    x.eval.grid <- obj$x.eval.grid
  }
  
  if ( mode == "df" ){
    tf.fit <- trendfilter(data$x, 
                          data$y, 
                          data$weights, 
                          k = obj$k,
                          lambda = obj$lambda, 
                          control = trendfilter.control.list(max_iter = obj$max_iter,
                                                             obj_tol = obj$obj_tol
                                                             )
                          )
    
    i.min <- which.min( abs(tf.fit$df - obj$df.min) )
    lambda.min <- obj$lambda[i.min]
    
    if ( obj$df[i.min] <= 2 ){
      return(NULL)
    }
    
  }else{
    tf.fit <- trendfilter(data$x, 
                          data$y, 
                          data$weights,
                          k = obj$k,
                          lambda = obj$lambda.min, 
                          control = glmgen::trendfilter.control.list(max_iter = obj$max_iter, 
                                                                     obj_tol = obj$obj_tol
                                                                     )
                          )
    
    lambda.min <- obj$lambda.min
  }

  tf.estimate <- as.numeric(predict(tf.fit, x.new = x.eval.grid, lambda = lambda.min))
  
  return(tf.estimate)
}


#' @importFrom dplyr slice_sample
nonparametric.resampler <- function(data){
  resampled.data <- slice_sample(data, n = nrow(data), replace = TRUE)
  return(resampled.data)
}


#' @importFrom stats rnorm
parametric.sampler <- function(data){
  boot.sample <- data$tf.estimate + rnorm(nrow(data), sd = 1 / sqrt(data$weights))
  return(data.frame(x = data$x, y = boot.sample, weights = data$weights))
}


wild.sampler <- function(data){
  wild.boot.residuals <- data$residuals * sample(x = c((1+sqrt(5))/2, 1-sqrt(5)/2), 
                                                 size = nrow(data), 
                                                 replace = TRUE, 
                                                 prob = c((1+sqrt(5))/(2*sqrt(5)), (sqrt(5)-1)/(2*sqrt(5)))
                                                 )
  wild.boot.sample <- data$tf.estimate + wild.boot.residuals
  return(data.frame(x = data$x, y = wild.boot.sample, weights = data$weights))
}
