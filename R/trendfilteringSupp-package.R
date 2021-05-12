#' Optimal one-dimensional data analysis with trend filtering
#'
#' @description This package serves as a software supplement to 
#' Politsch et al. (2020a) and Politsch et al. (2020b). We build on 
#' the glmgen package by providing additional useful statistical tools 
#' for carrying out data analysis with trend filtering 
#' (Tibshirani 2014). It contains user-friendly functionality for 
#' optimizing a trend filtering estimator by cross validation or 
#' minimizing Stein's unbiased risk estimate, as well as various 
#' bootstrap algorithms for producing variability bands to quantify 
#' the uncertainty in the estimator.
#' @name trendfilteringSupp-package
#' @docType package
#' @author Collin A. Politsch \cr \cr 
#' \strong{Maintainer}: Collin A. Politsch <collinpolitsch@@gmail.com>
#' @keywords package
#' @references 
#' \enumerate{
#' \item{Politsch et al. (2020a). Trend filtering – I. A modern statistical tool
#' for time-domain astronomy and astronomical spectroscopy. \emph{Monthly 
#' Notices of the Royal Astronomical Society}, 492(3), p. 4005-4018. 
#' \href{https://academic.oup.com/mnras/article/492/3/4005/5704413}{[Link]}} \cr
#' \item{Politsch et al. (2020b). Trend Filtering – II. Denoising astronomical 
#' signals with varying degrees of smoothness. \emph{Monthly Notices of the 
#' Royal Astronomical Society}, 492(3), p. 4019-4032.
#' \href{https://academic.oup.com/mnras/article/492/3/4019/5704414}{[Link]}} \cr
#' \item{Tibshirani (2014). Adaptive piecewise polynomial estimation via trend 
#' filtering. \emph{The Annals of Statistics}. 42(1), p. 285-323.
#' \href{https://projecteuclid.org/euclid.aos/1395234979}{[Link]}} \cr
#' \item{Ramdas and Tibshirani (2016). Fast and Flexible ADMM Algorithms 
#' for Trend Filtering. \emph{Journal of Computational and Graphical 
#' Statistics}, 25(3), p. 839-858.
#' \href{https://amstat.tandfonline.com/doi/abs/10.1080/10618600.2015.1054033#.XfJpNpNKju0}{[Link]}} \cr
#' \item{Arnold, Sadhanala, and Tibshirani (2014). Fast algorithms for 
#' generalized lasso problems. R package \emph{glmgen}. Version 0.0.3. 
#' \href{https://github.com/glmgen/glmgen}{[Link]}} \cr
#' \item{Tibshirani and Taylor (2012). Degrees of freedom in lasso problems.
#' \emph{The Annals of Statistics}, 40(2), p. 1198-1232.
#' \href{https://projecteuclid.org/journals/annals-of-statistics/volume-40/issue-2/Degrees-of-freedom-in-lasso-problems/10.1214/12-AOS1003.full}{[Link]}} \cr
#' \item{Efron and Tibshirani (1986). Bootstrap Methods for Standard Errors, 
#' Confidence Intervals, and Other Measures of Statistical Accuracy. Statistical
#' Science, 1(1), p. 54-75.
#' \href{https://projecteuclid.org/journals/statistical-science/volume-1/issue-1/Bootstrap-Methods-for-Standard-Errors-Confidence-Intervals-and-Other-Measures/10.1214/ss/1177013815.full}{[Link]}} \cr
#' \item{Wu (1986). Jackknife, Bootstrap and Other Resampling Methods in 
#' Regression Analysis. \emph{The Annals of Statistics}, 14(4), 1261-1295.
#' \href{https://projecteuclid.org/journals/annals-of-statistics/volume-14/issue-4/Jackknife-Bootstrap-and-Other-Resampling-Methods-in-Regression-Analysis/10.1214/aos/1176350142.full}{[Link]}} \cr
#' \item{Efron (1979). Bootstrap Methods: Another Look at the Jackknife.
#' \emph{The Annals of Statistics}, 7(1), p. 1-26.
#' \href{https://projecteuclid.org/journals/annals-of-statistics/volume-7/issue-1/Bootstrap-Methods-Another-Look-at-the-Jackknife/10.1214/aos/1176344552.full}{[Link]}} \cr
#' \item{Tibshirani and Wasserman (2015). Stein’s Unbiased Risk Estimate.
#' \emph{36-702: Statistical Machine Learning course notes} (Carnegie Mellon).
#' \href{http://www.stat.cmu.edu/~larry/=sml/stein.pdf}{[Link]}} \cr
#' \item{Efron (2014). The Estimation of Prediction Error: Covariance Penalties 
#' and Cross-Validation. \emph{Journal of the American Statistical Association}.
#' 99(467), p. 619-632.
#' \href{https://www.tandfonline.com/doi/abs/10.1198/016214504000000692}{[Link]}} \cr
#' \item{Stein (1981). Estimation of the Mean of a Multivariate Normal 
#' Distribution. \emph{The Annals of Statistics}. 9(6), p. 1135-1151.
#' \href{https://projecteuclid.org/journals/annals-of-statistics/volume-9/issue-6/Estimation-of-the-Mean-of-a-Multivariate-Normal-Distribution/10.1214/aos/1176345632.full}{[Link]}} \cr
#' \item{Hastie, Tibshirani, and Friedman (2009). The Elements of Statistical 
#' Learning: Data Mining, Inference, and Prediction. 2nd edition. Springer 
#' Series in Statistics. \href{https://web.stanford.edu/~hastie/ElemStatLearn/printings/ESLII_print12_toc.pdf}{
#' [Online print #12]}} \cr
#' }
NULL

