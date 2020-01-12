\name{SURE.trendfilter}
\alias{SURE.trendfilter}
\title{SURE.trendfilter}
\usage{
SURE.trendfilter(x = NULL, y, sigma = NULL, lambda = NULL, k = 2, max_iter = 250, obj_tol = 1e-06)
}
\arguments{
\item{x}{The observed inputs.}
\item{y}{The observed outputs.}
\item{sigma}{A vector of measurement standard errors of the observed y's.}
\item{lambda}{A vector of trend filtering hyperparameter values to be tested.}
\item{k}{The degree of the trend filtering estimator. Defaults to \code{k=2} (quadratic trend filtering).}
\item{max_iter}{Maximum iterations for ADMM optimization. Consider increasing if estimate does not appear to have fully converged.}
\item{obj_tol}{The tolerance used in the stopping criterion; when the relative change in objective values is less than this value, the algorithm terminates.}
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
install.packages("devtools")
devtools::install_github("statsmaths/glmgen", subdir="R_pkg/glmgen")
  
# Quasar spectrum example
data(quasar)
  
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
lambda.grid <- exp(seq(-10,7,length=250))
SURE.out <- SURE.trendfilter(log.wavelength.scaled, flux, wts, lambda.grid)
lambda.opt <- SURE.out$lambda[which.min(SURE.out$SURE.loss)]

# Fit optimized model
fit <- glmgen::trendfilter(log.wavelength.scaled, flux, wts, k = 2, lambda = lambda.opt)

# Plot results
par(mfrow=c(2,1))
plot(log(lambda.grid), SURE.out$SURE.loss, xlab = "log(lambda)", ylab = "SURE")
abline(v = log(lambda.opt))
wavelength <- 10 ^ (log.wavelength.scaled / 1000)
plot(wavelength, flux, type = "l")
lines(wavelength, fit$beta, col = "orange", lwd = 2.5)
}
