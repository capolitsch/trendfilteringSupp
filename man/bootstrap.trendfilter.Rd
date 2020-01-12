\name{bootstrap.trendfilter}
\alias{bootstrap.trendfilter}
\title{bootstrap.trendfilter}
\usage{
bootstrap.trendfilter(
  x,
  y,
  lambda.opt,
  sigma = NULL,
  B = 1000,
  x.eval.grid = x,
  k = 2,
  bootstrap.method = "nonparametric",
  alpha = 0.05,
  full.ensemble = F,
  max_iter = 250,
  obj_tol = 1e-06,
  mc.cores = 1
)
}
\arguments{
\item{x}{The observed inputs.}
\item{y}{The observed outputs.}
\item{lambda.opt}{The optimally-tuned trend filtering hyperparameter, e.g. by cross validation on minimizing SURE (see \code{\link{SURE.trendfilter}}).}
\item{sigma}{A vector of measurement standard errors of the observed y's.}
\item{B}{The number of boostrap samples used to compute the confidence bands.}
\item{x.eval.grid}{Input evaluation grid. Defaults to the observed inputs.}
\item{k}{The degree of the trend filtering estimator. Defaults to \code{k=2} (quadratic trend filtering).}
\item{bootstrap.method}{Bootstrap method to be implemented. See description for suggested use. Defaults to "nonparametric".}
\item{alpha}{Specifies the width of the \code{1-alpha} pointwise confidence bands.}
\item{full.ensemble}{Return the full bootstrap ensemble as an (n x B) matrix. Defaults to FALSE.}
\item{max_iter}{Maximum iterations for ADMM optimization. Consider increasing if estimate does not appear to have fully converged.}
\item{obj_tol}{The tolerance used in the stopping criterion; when the relative change in objective values is less than this value, the algorithm terminates.}
\item{mc.cores}{The number of cores to use.}
}

\value{
A list with the following elements:
\item{bootstrap.lower.perc.intervals}{Vector of lower bounds for the \code{1-alpha} pointwise confidence band.}
\item{bootstrap.upper.perc.intervals}{Vector of upper bounds for the \code{1-alpha} pointwise confidence band.}
\item{tf.boot.ensemble}{The full bootstrap ensemble as an (n x B) matrix.}
}
\description{
bootstrap.trendfilter takes data of the form (x_1,y_1), ...,(x_n,y_n) or (x_1,y_1,sigma_1), ...,(x_n,y_n,sigma_n) and implements one of three bootstrap algorithms. The bootstrap method should generally be chosen according to the following criteria: 
S1. The inputs are irregularly sampled --> \code{bootstrap.method = "nonparametric"};
S2. The inputs are regularly sampled and the noise distribution is known --> \code{bootstrap.method = "parametric"};
S3. The inputs are regularly sampled and the noise distribution is unknown --> \code{bootstrap.method = "wild"}.
}
\author{
Collin Politsch <collinpolitsch@gmail.com>
}
\seealso{
  \code{\link{SURE.trendfilter}}
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
wavelength <- 10 ^ (log.wavelength.scaled / 1000)
plot(wavelength, flux, type = "l")
lines(wavelength, fit$beta, col = "orange", lwd = 2.5)

boot.out <- bootstrap.trendfilter(log.wavelength.scaled, flux, lambda.opt, sigma = sqrt(1/wts),
                                  bootstrap.method = "parametric")
lines(wavelength, boot.out$bootstrap.lower.perc.intervals, col = "orange", lty = 2, lwd = 2)
lines(wavelength, boot.out$bootstrap.upper.perc.intervals, col = "orange", lty = 2, lwd = 2)
legend(x = "topleft", lty = c(1,2), col = "orange", lwd = 2, 
       legend = c("Trend filtering estimate", "95 percent variability band"))
}
