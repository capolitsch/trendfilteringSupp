library(devtools)
install_github("statsmaths/glmgen", subdir="R_pkg/glmgen")

SURE.trendfilter <- function(x = NULL, y, wts = NULL, gamma = NULL, deg = 2){
  if ( is.null(x) ){
    x <- rep(1, length(y))
  }
  if ( is.null(wts) ){
    stop("Weights must be provided in order to compute SURE.")
  }
  if ( is.null(gamma) ){
    stop("lambda must be specified.")
  }
  out <- glmgen::trendfilter(x = x, y = y, weights = wts, k = deg, lambda = gamma)
  if ( length(gamma) == 1 ){
    SURE.loss <- mean( (out$beta - y)^2 ) + (2 * mean(1/wts) / length(x)) * out$df
  }
  if ( length(gamma) > 1 ){
    SURE.loss <- colMeans( (out$beta - y)^2 ) + (2 * mean(1/wts) / length(x)) * out$df
  }
  return(list(gamma = gamma, SURE.loss = as.numeric(SURE.loss)))
}

# Quasar spectrum example
install.packages("FITSio")
library(FITSio)
quasar.spec <- readFITS("spec-4055-55359-0010.fits", hdu = 1)

# SDSS spectra are equally spaced in log10 wavelength space with a separation of 10e-4
# Reading in a spectrum file and retrieving the piece of the spectrum in the Lyman-alpha
# forest region
log.wavelength.scaled <- quasar.spec$col[[2]] * 1000
flux <- quasar.spec$col[[1]]
wts <- quasar.spec$col[[3]]
lya.rest.wavelength <- 1215.67
inds <- which(( 10 ^ (log.wavelength.scaled / 1000) ) / (2.953 + 1) < lya.rest.wavelength + 40)
log.wavelength.scaled <- log.wavelength.scaled[inds]
flux <- flux[inds]
wts <- wts[inds]

# Compute SURE loss curve and optimal gamma
gamma.grid <- exp(seq(-14,4,length=150))
SURE.out <- SURE.trendfilter(log.wavelength.scaled, flux, wts, gamma.grid)
gamma.opt <- SURE.out$gamma[which.min(SURE.out$SURE.loss)]

# Fit optimized model
fit <- glmgen::trendfilter(log.wavelength.scaled, flux, wts, k = 2, lambda = gamma.opt)

# Plot results
wavelength <- 10 ^ (log.wavelength.scaled / 1000)
plot(wavelength, flux, type = "l")
lines(wavelength, fit$beta, col = "orange", lwd = 2.5)
