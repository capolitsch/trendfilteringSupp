\name{SURE.trendfilter}
\alias{SURE.trendfilter}
\title{SURE.trendfilter}
\usage{
SURE.trendfilter(x = NULL, y, sigma = NULL, lambda = NULL, k = 2)
}
\arguments{
\item{x}{The observed inputs.}
\item{y}{The observed outputs.}
\item{sigma}{A vector of measurement standard errors of the observed y's.}
\item{lambda}{A vector of trend filtering hyperparameter values to be tested.}
\item{k}{The degree of the trend filtering estimator. Defaults to \code{k=2} (quadratic trend filtering).}
}

\value{
A list with the following elements:
\item{lambda}{Vector of hyperparameter values tested.}
\item{SURE.loss}{Vector of estimated SURE errors for each hyperparameter value.}
}
\description{
\code{SURE.trendfilter} computes the Stein's unbiased risk estimate of fixed-input mean-squared prediction error on a grid of hyperparameter values.
}
\author{
Collin Politsch <collinpolitsch@gmail.com>
}
\seealso{
  \code{\link{bootstrap.trendfilter}}
}
\examples{
# Quasar spectrum example
install.packages("FITSio")
install.packages("devtools")
library(FITSio)
library(devtools)
install_github("statsmaths/glmgen", subdir="R_pkg/glmgen")
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

# Compute SURE loss curve and optimal lambda
lambda.grid <- exp(seq(-14,4,length=150))
SURE.out <- SURE.trendfilter(log.wavelength.scaled, flux, wts, lambda.grid)
lambda.opt <- SURE.out$lambda[which.min(SURE.out$SURE.loss)]

# Fit optimized model
fit <- glmgen::trendfilter(log.wavelength.scaled, flux, wts, k = 2, lambda = lambda.opt)

# Plot results
wavelength <- 10 ^ (log.wavelength.scaled / 1000)
plot(wavelength, flux, type = "l")
lines(wavelength, fit$beta, col = "orange", lwd = 2.5)
}
