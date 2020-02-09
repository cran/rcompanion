#' @title        Quantiles and confidence intervals 
#'
#' @description  Calculates an estimate for a quantile 
#'               and confidence intervals for
#'               a vector of discrete or continuous
#'               values
#'               
#' @param x      The vector of observations.
#'               Can be an ordered factor as long as \code{type}
#'               is 1 or 3.
#' @param tau    The quantile to use,
#'               e.g. 0.5 for median, 0.25 for 25th percentile.
#' @param level  The confidence interval to use, 
#'               e.g. 0.95 for 95 percent confidence interval. 
#' @param method If \code{"binomial"}, uses the binomial distribution 
#'               the confidence limits.
#'               If \code{"normal"}, uses the normal approximation to the
#'               binomial distribution.
#' @param type   The \code{type} value passed to the \code{quantile} function.
#' @param digits The number of significant figures to use in output.
#' @param ...    Other arguments, ignored.
#' 
#' @details      Conover recommends the \code{"binomial"} method for sample
#'               sizes less than or equal to 20.
#'               With the current implementation, 
#'               this method can be used also for
#'               larger sample sizes.
#'               
#' @author       Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references   \url{http://rcompanion.org/handbook/E_04.html},
#'               Conover, W.J., Practical Nonparametric Statistics, 3rd.
#' @seealso      \code{\link{groupwisePercentile}}, \code{\link{groupwiseMedian}}
#' @concept      percentile quantile confidence interval
#' @return       A data frame of summary statistics, quantile estimate, 
#'               and confidence limits.
#'          
#' @examples
#' ### From Conover, Practical Nonparametric Statistics, 3rd
#' Hours = c(46.9, 47.2, 49.1, 56.5, 56.8, 59.2, 59.9, 63.2,
#'           63.3, 63.4, 63.7, 64.1, 67.1, 67.7, 73.3, 78.5)
#' quantileCI(Hours)
#' 
#' ### Example with ordered factor
#' set.seed(12345)
#' Pool = factor(c("smallest", "small", "medium", "large", "largest"),
#'              ordered=TRUE, 
#'              levels=c("smallest", "small", "medium", "large", "largest"))
#' Sample = sample(Pool, 24, replace=TRUE)
#' quantileCI(Sample)
#'                  
#' @importFrom stats quantile qbinom qnorm pbinom
#' 
#' @export

quantileCI = function(x, tau=0.50, level=0.95, method="binomial",
                      type=3, digits=3, ...){
  n       = length(x)
  q       = tau
  Ordered = sort(x)
  
  if(method=="binomial"){
     lwr     = qbinom((1-level)/2, n, tau)
     upr     = qbinom(1-(1-level)/2, n, tau)
     LWR     = Ordered[lwr]
     UPR     = Ordered[upr]
  }
  
  if(method=="normal"){
     p       = qnorm((1-level)/2, lower.tail=FALSE)
     lwr     =   floor(n*q - p*sqrt(n*q*(1-q)))
     upr     = ceiling(n*q + p*sqrt(n*q*(1-q)))
     LWR     = Ordered[lwr]
     UPR     = Ordered[upr]
  }
  
  QUANT      = quantile(Ordered, q, type=type)
  LEVEL1     = pbinom(lwr,n,tau)
  LEVEL2     = pbinom(upr,n,tau)
  ACTUAL     = LEVEL2 - LEVEL1

  if(is.numeric(x)){
     Z          = data.frame(tau           = tau, 
                             n             = n, 
                             Quantile      = signif(QUANT, digits=digits), 
                             Nominal.level = level,
                             Actual.level  = signif(ACTUAL, digits=digits),
                             Lower.ci      = signif(LWR, digits=digits), 
                             Upper.ci      = signif(UPR, digits=digits))
  }
  
  if(is.factor(x)){
     Z          = data.frame(tau           = tau, 
                             n             = n, 
                             Quantile      = QUANT, 
                             Nominal.level = level,
                             Actual.level  = signif(ACTUAL, digits=digits),
                             Lower.ci      = LWR, 
                             Upper.ci      = UPR)
     }
  row.names(Z) = ""
  return(Z)
}