# Optimal one-dimensional data analysis with trend filtering

```
devtools::install_github("capolitsch/trendfilteringSupp")
library(trendfilteringSupp)

?SURE.trendfilter
?cv.trendfilter
?bootstrap.trendfilter
```

This package serves as a software supplement to [Politsch et al. (2020a)](https://academic.oup.com/mnras/article/492/3/4005/5704413) 
and [Politsch et al. (2020b)](https://academic.oup.com/mnras/article/492/3/4019/5704414).
We provide a variety of statistical tools for one-dimensional data analyses 
with trend filtering [(Tibshirani 2014)](https://projecteuclid.org/euclid.aos/1395234979). 
This package contains user-friendly functionality for optimizing a trend 
filtering estimator by cross validation or Stein's unbiased risk estimate and 
various bootstrap algorithms for producing variability bands to quantify the 
uncertainty in the optimized trend filtering estimate.


## References:

Politsch et al. Trend Filtering – I. A modern statistical tool for time-domain astronomy and Astronomical Spectroscopy. 
*Monthly Notices of the Royal Astronomical Society*, 492(3), p. 4005-4018, 2020. [[Link](https://academic.oup.com/mnras/article/492/3/4005/5704413)]

Politsch et al. Trend Filtering – II. Denoising astronomical signals with varying degrees of smoothness. 
*Monthly Notices of the Royal Astronomical Society*, 492(3), p. 4019-4032, 2020. [[Link](https://academic.oup.com/mnras/article/492/3/4019/5704414)]

R. J. Tibshirani. Adaptive piecewise polynomial estimation via trend filtering. 
*The Annals of Statistics*. 42 (2014), no. 1, 285-323. [[Link](https://projecteuclid.org/euclid.aos/1395234979)]

A. Ramdas & R. J. Tibshirani. Fast and Flexible ADMM Algorithms for Trend Filtering.
*Journal of Computational and Graphical Statistics*, 25:3 (2016), 839-858. [[Link](https://amstat.tandfonline.com/doi/abs/10.1080/10618600.2015.1054033#.XfJpNpNKju0)]

T. B. Arnold, V. Sadhanala, and R. J. Tibshirani. Fast algorithms for generalized lasso problems. *glmgen* R package.
version 0.0.3 (2014). [[Link](https://github.com/glmgen/glmgen)]


## Quasar Lyman-alpha forest example 

```
# Load Lyman-alpha forest spectral observations of an SDSS quasar at redshift 
# z ~ 2.953. SDSS spectra are equally spaced in log10 wavelength space, 
# aside from some instances of masked pixels.
 
 data(quasar_spec)
 data(plotting_utilities)
 
 
 # We are interested in denoising the observed brightness of the quasar 
 # (measured as a 'flux' quantity) over the observed wavelength range. Since 
 # the logarithmic wavelengths are gridded, we optimize the trend filtering 
 # hyperparameter by minimizing the SURE estimate of fixed-input squared 
 # prediction error. For smoothness, we use quadratic trend filtering, i.e. 
 # the default k=2. 
 
 SURE.obj <- SURE.trendfilter(x = log10.wavelength, 
                              y = flux, 
                              weights = weights)
 
 
 # Extract the SURE error curve and optimized trend filtering estimate from 
 # the SURE.trendfilter output
 
 log.lambda <- log(SURE.obj$lambda)
 error <- SURE.obj$error
 log.lambda.min <- log(SURE.obj$lambda.min)
 
 log10.wavelength.eval <- SURE.obj$x.eval
 tf.estimate <- SURE.obj$tf.estimate
 
 
 # Run a parametric bootstrap on the optimized trend filtering estimator to 
 # obtain uncertainty bands
 
 boot.out <- bootstrap.trendfilter(obj = SURE.obj, 
                                   bootstrap.algorithm = "parametric")
 
 
 # Transform the inputs to wavelength space (in Angstroms)
 
 wavelength <- 10 ^ (log10.wavelength)
 wavelength.eval <- 10 ^ (log10.wavelength.eval)
 
 
 # Plot the results

 par(mfrow = c(2,1), mar = c(5,4,2.5,1) + 0.1)
 plot(x = log.lambda, y = error, main = "SURE error curve", 
      xlab = "log(lambda)", ylab = "SURE error")
 abline(v = log.lambda.min, lty = 2, col = "blue3")
 text(x = log.lambda.min, y = par("usr")[4], 
      labels = "optimal hyperparameter", pos = 1, col = "blue3")
 
 plot(x = wavelength, y = flux, type = "l", 
      main = "Quasar Lyman-alpha forest", 
      xlab = "Observed wavelength (Angstroms)", ylab = "Flux")
 lines(wavelength.eval, tf.estimate, col = "orange", lwd = 1.5)
 polygon(c(wavelength.eval, rev(wavelength.eval)), 
         c(boot.out$bootstrap.lower.band, 
         rev(boot.out$bootstrap.upper.band)),
         col = transparency("orange", 90), border = NA)
 lines(wavelength.eval, boot.out$bootstrap.lower.band, 
       col = "orange", lwd = 0.5)
 lines(wavelength.eval, boot.out$bootstrap.upper.band, 
       col = "orange", lwd = 0.5)
 legend(x = "topleft", lwd = c(1,2,8), lty = 1, cex = 0.75,
        col = c("black","orange", transparency("orange", 90)), 
        legend = c("Noisy quasar spectrum",
                   "Trend filtering estimate",
                   "95% variability band"))
```
