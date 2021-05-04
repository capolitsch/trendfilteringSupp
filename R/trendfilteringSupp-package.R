#' Supplementary tools for data analysis with trend filtering
#'
#' @description This package builds on \code{glmgen} by providing additional 
#'              useful statistical tools for carrying out data analysis with 
#'              trend filtering [Tibshirani (2014)]. It contains functionality 
#'              for easily optimizing a trend filtering estimator by cross 
#'              validation or minimizing Stein's unbiased risk estimate and 
#'              various bootstrap algorithms for producing variability bands to 
#'              quantify the uncertainty in the estimator.
#' @name trendfilteringSupp-package
#' @docType package
#' @author Collin A. Politsch \cr \cr 
#' \strong{Maintainer}: Collin A. Politsch <collinpolitsch@@gmail.com>
#' @import glmgen, dplyr
#' @keywords package
#' @seealso Refer to:
#' \enumerate{
#' \item{\href{https://academic.oup.com/mnras/article/492/3/4005/5704413}{Trend filtering – I. A modern statistical tool for time-domain astronomy and astronomical spectroscopy}}
#' \item{\href{https://academic.oup.com/mnras/article/492/3/4019/5704414}{Trend filtering – II. Denoising astronomical signals with varying degrees of smoothness}}
#' }
NULL

