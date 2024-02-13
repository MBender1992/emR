#' Calculate pvalues from confidence intervals
#'
#' The calculation method is based on this publication: https://www.bmj.com/content/343/bmj.d2304.
#' @param est Estimate (e.g. difference or ratio)
#' @param l lower limit of confidence interval
#' @param u upper limit of confidence interval
#' @param log.trans logical indicating whether estimates should be log transformed. If TRUE values will be log transformed before calculation. Default is TRUE.
#' @param conf.level confidence level
#' @export

ci_pval <- function(est, l, u, log.trans = TRUE, conf.level = 0.95){
  if(log.trans == TRUE){
    est <- log(est)
    l   <- log(l)
    u   <- log(u)
  }
  quantile <- (1-conf.level)/2
  SE <- (u - l)/(2*stats::qnorm(quantile,lower.tail=FALSE))
  z <- abs(est/SE)
  p <- exp(-0.717*z - 0.416*z^2)
  return(p)
}


