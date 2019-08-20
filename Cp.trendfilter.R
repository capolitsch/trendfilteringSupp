library(devtools)
install_github("statsmaths/glmgen", subdir="R_pkg/glmgen")

Cp.trendfilter <- function(x = NULL, y, wts = NULL, gamma = NULL, deg = 2){
  if ( is.null(x) ){
    x <- rep(1, length(y))
  }
  if ( is.null(wts) ){
    stop("Weights must be provided in order to compute Mallows' Cp.")
  }
  if ( is.null(gamma) ){
    stop("lambda must be specified.")
  }
  out <- glmgen::trendfilter(x = x, y = y, weights = wts, k = deg, lambda = gamma)
  if ( length(gamma) == 1 ){
    Cp.loss <- mean( (out$beta - y)^2 ) + (2 * mean(1/wts) / length(x)) * out$df
  }
  if ( length(gamma) > 1 ){
    Cp.loss <- colMeans( (out$beta - y)^2 ) + (2 * mean(1/wts) / length(x)) * out$df
  }
  return(list(gamma = gamma, Cp.loss = as.numeric(Cp.loss)))
}

# Galaxy spectrum example
install.packages("FITSio")
library(FITSio)
galaxy.spec <- readFITS("spec-4055-55359-0018.fits", hdu = 1)

log.wavelength.scaled <- galaxy.spec$col[[2]] * 1000
flux <- galaxy.spec$col[[1]]
wts <- galaxy.spec$col[[3]]

gamma.grid <- exp(seq(-10,6,length=100))
Cp.out <- Cp.trendfilter(log.wavelength.scaled, flux, wts, gamma.grid)
gamma.opt <- Cp.out$gamma[which.min(Cp.out$Cp.loss)]
fit <- glmgen::trendfilter(log.wavelength.scaled, flux, wts, k = 2, lambda = gamma.opt)

wavelength <- 10 ^ (log.wavelength.scaled / 1000)
plot(wavelength, flux, type = "l")
lines(wavelength, fit$beta, col = "orange", lwd = 2)
