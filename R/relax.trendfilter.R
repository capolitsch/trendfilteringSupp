#' Relaxed trend filtering (under construction. do not use.)
#'
#' @description \loadmathjax{} \code{relax.trendfilter} 
#' @param gamma mixing hyperparameter
#' @return 
#' @details This should be a very detailed description...
#' @export relax.trendfilter
#' @author Collin A. Politsch, \email{collinpolitsch@@gmail.com}
#' @seealso \code{\link{trendfilter}}, \code{\link{SURE.trendfilter}}, 
#' \code{\link{cv.trendfilter}}, \code{\link{bootstrap.trendfilter}}
#' @references \enumerate{
#' \item{Politsch et al. (2020a). Trend filtering – I. A modern 
#' statistical tool for time-domain astronomy and astronomical spectroscopy. 
#' \emph{Monthly Notices of the Royal Astronomical Society}, 492(3), 
#' p. 4005-4018.
#' \href{https://academic.oup.com/mnras/article/492/3/4005/5704413}{[Link]}} \cr
#' 
#' \item{Politsch et al. (2020b). Trend Filtering – II. Denoising 
#' astronomical signals with varying degrees of smoothness. \emph{Monthly 
#' Notices of the Royal Astronomical Society}, 492(3), p. 4019-4032.
#' \href{https://academic.oup.com/mnras/article/492/3/4019/5704414}{[Link]}} \cr
#' 
#' \item{Meinshausen (2007). Relaxed Lasso. \emph{Computational Statistics & 
#' Data Analysis}, 52(1), p. 374-393.
#' \href{https://www.sciencedirect.com/science/article/abs/pii/S0167947306004956}{[Link]}} \cr
#' }

relax.trendfilter <- function(obj, gamma
                              )
{
  mc.cores <- max(c(parallel::detectCores() - 2), 1)
  
}