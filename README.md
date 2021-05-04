# Statistical analysis with trend filtering
This package serves as a supplement to [Politsch et al. (2020a)](https://academic.oup.com/mnras/article/492/3/4005/5704413) 
and [Politsch et al. (2020b)](https://academic.oup.com/mnras/article/492/3/4019/5704414).
The package builds on [*glmgen*](https://github.com/glmgen/glmgen) by providing 
additional useful statistical tools for carrying out data analysis with trend 
filtering [Tibshirani (2014)](https://projecteuclid.org/euclid.aos/1395234979). 
It contains functionality for easily optimizing a trend filtering estimator by 
cross validation or minimizing Stein's unbiased risk estimate and various 
bootstrap algorithms for producing variability bands to quantify the uncertainty 
in the estimator.

```
devtools::install_github("capolitsch/trendfilteringSupp")
library(trendfilteringSupp)

?SURE.trendfilter
?bootstrap.trendfilter
```

# Quasar Lyman-alpha forest example 

SDSS spectra are equally spaced in log base-10 wavelength space with a 
separation of 10e-4 log-Angstroms. Given the default trend filtering 
optimization parameters, it is safer to scale up the inputs in such a 
scenario. For example, here we scale to unit spacing.

```
# Read in an SDSS spectrum of a quasar at redshift z = 2.953 and extract the 
# Lyman-alpha forest.

data(quasar_spec)
lya.rest <- 1215.67
quasar.redshift <- 2.953

log.wavelength.scaled <- quasar_spec$col[[2]] * 1000
flux <- quasar_spec$col[[1]]
weights <- quasar_spec$col[[3]]

inds <- which((10^(quasar_spec$col[[2]]))/(quasar.redshift + 1) < lya.rest)
x <- log.wavelength.scaled[inds]
y <- flux[inds]
weights <- weights[inds]
# Run the SURE optimization for a quadratic trend filtering estimator, i.e. 
# k = 2 (recommended)

set.seed(1)
lambda.grid <- exp(seq(-10, 5, length = 200))
SURE.obj <- SURE.trendfilter(x = x, 
                            y = y, 
                            weights = weights, 
                            k = 2,
                            lambda = lambda.grid
                            )
lambda.min <- SURE.obj$lambda.min


# Fit the optimized trend filtering model and get the estimates on an fine
# equally-spaced input grid

model <- trendfilter(x = x,
                    y = y, 
                    weights = weights,
                    k = 2, 
                    lambda = lambda.min
                    )
                    
x.eval <- seq(min(x), max(x), length = 1500)
tf.estimate <- predict(model, x.new = x.eval)


# Plot the results
par(mfrow = c(2,1), mar = c(5,4,2.5,1) + 0.1)

plot(log(lambda.grid), SURE.obj$error,
    main = "SURE error curve", 
    xlab = "log(lambda)", ylab = "SURE error")
abline(v = log(lambda.min), col = "blue3", lty = 2)
text(x = log(lambda.min), y = par("usr")[4], 
    labels = "optimal hyperparameter", pos = 1, col = "blue3")
    
    
# Transform back to wavelength space
wavelength <- 10 ^ (x / 1000)
wavelength.eval <- 10 ^ (x.eval / 1000)

plot(wavelength, y, type = "l", 
    main = "Quasar Lyman-alpha forest", 
    xlab = "Observed wavelength (angstroms)", ylab = "flux")
lines(wavelength.eval, tf.estimate, col = "orange", lwd = 2.5)
legend(x = "topleft", lwd = c(1,2), lty = 1, 
      col = c("black","orange"), 
      legend = c("Noisy quasar spectrum",
                 "Trend filtering estimate"
                 )
      )
```

References:

C. A. Politsch et al. Trend Filtering - I. A modern statistical tool for time-domain astronomy 
and Astronomical Spectroscopy. Monthly Notices of the Royal Astronomical Society, 492(3), p. 4005-4018, 2020. 
https://academic.oup.com/mnras/article/492/3/4005/5704413

C. A. Politsch et al. Trend Filtering - II. Denoising astronomical signals with varying degrees of smoothness. Monthly Notices of the Royal Astronomical Society, 492(3), p. 4019-4032, 2020. 
https://academic.oup.com/mnras/article/492/3/4019/5704414

R. J. Tibshirani. Adaptive piecewise polynomial estimation via trend filtering. 
The Annals of Statistics. 42 (2014), no. 1, 285--323. doi:10.1214/13-AOS1189. 
https://projecteuclid.org/euclid.aos/1395234979

A. Ramdas & R. J. Tibshirani. Fast and Flexible ADMM Algorithms for Trend Filtering.
Journal of Computational and Graphical Statistics, 25:3 (2016), 839-858, DOI: 10.1080/10618600.2015.1054033.
https://amstat.tandfonline.com/doi/abs/10.1080/10618600.2015.1054033#.XfJpNpNKju0

T. B. Arnold, V. Sadhanala, and R. J. Tibshirani. Fast algorithms for generalized lasso problems.
https://github.com/glmgen. Version 0.0.3 (2014)
