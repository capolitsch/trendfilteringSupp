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


## Key references:

C. A. Politsch et al. (2020a). Trend Filtering – I. A modern statistical tool for time-domain astronomy and Astronomical Spectroscopy. 
*Monthly Notices of the Royal Astronomical Society*, 492(3), p. 4005-4018. [[Link](https://academic.oup.com/mnras/article/492/3/4005/5704413)]

C. A. Politsch et al. (2020b). Trend Filtering – II. Denoising astronomical signals with varying degrees of smoothness. 
*Monthly Notices of the Royal Astronomical Society*, 492(3), p. 4019-4032. [[Link](https://academic.oup.com/mnras/article/492/3/4019/5704414)]

R. J. Tibshirani (2014). Adaptive piecewise polynomial estimation via trend filtering. 
*The Annals of Statistics*, 42(1), p. 285-323. [[Link](https://projecteuclid.org/euclid.aos/1395234979)]

A. Ramdas & R. J. Tibshirani (2016). Fast and Flexible ADMM Algorithms for Trend Filtering.
*Journal of Computational and Graphical Statistics*, 25(3), p. 839-858.
[[Link](https://amstat.tandfonline.com/doi/abs/10.1080/10618600.2015.1054033#.XfJpNpNKju0)]

T. B. Arnold, V. Sadhanala, and R. J. Tibshirani (2014). Fast algorithms for generalized lasso problems. *glmgen* R package,
version 0.0.3. [[Link](https://github.com/glmgen/glmgen)]
