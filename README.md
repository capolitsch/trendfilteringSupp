# Trend filtering
Supplementary tools for data analysis with trend filtering (wrapper of glmgen)

```
library(devtools)
install_github("statsmaths/glmgen", subdir="R_pkg/glmgen")
install_github("capolitsch/trendfilteringSupp")

library(glmgen)
library(trendfilteringSupp)
?SURE.trendfilter
?bootstrap.trendfilter
```

References:

C. A. Politsch et al. Trend Filtering - I. A Modern Statistical Tool for Time-Domain Astronomy 
and Astronomical Spectroscopy. Accepted for publication in Monthly Notices of the Royal Astronomical Society. https://arxiv.org/abs/1908.07151

C. A. Politsch et al. Trend Filtering - II. Denoising Astronomical Signals with Varying Degrees of Smoothness. Accepted for publication in Monthly Notices of the Royal Astronomical Society. https://arxiv.org/abs/2001.03552

Tibshirani, R. J. Adaptive piecewise polynomial estimation via trend filtering. 
The Annals of Statistics. 42 (2014), no. 1, 285--323. doi:10.1214/13-AOS1189. 
https://projecteuclid.org/euclid.aos/1395234979

A. Ramdas & R. J. Tibshirani. Fast and Flexible ADMM Algorithms for Trend Filtering.
Journal of Computational and Graphical Statistics, 25:3 (2016), 839-858, DOI: 10.1080/10618600.2015.1054033.
https://amstat.tandfonline.com/doi/abs/10.1080/10618600.2015.1054033#.XfJpNpNKju0

T. B. Arnold, V. Sadhanala, and R. J. Tibshirani. Fast algorithms for generalized lasso problems.
https://github.com/glmgen. Version 0.0.3 (2014)
