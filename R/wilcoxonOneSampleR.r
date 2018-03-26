#' @title r effect size for Wilcoxon one-sample signed-rank test
#' 
#' @description Calculates r effect size
#'              for a Wilcoxon one-sample signed-rank test.
#' 
#' @param x A vector of observations of an ordinal variable.
#' @param mu The value to compare \code{x} to, as in \code{wilcox.test}
#' @param digits The number of significant digits in the output.
#' @param ... Additional arguments passed to the \code{wilcoxsign_test} function.             
#' 
#' @details  A Z value is extracted from the \code{wilcoxsign_test} function in the
#'           coin package.  r  is calculated as Z divided by 
#'           square root of the number of observations.
#'           
#'           The calculated statistic is equivalent to the statistic returned
#'           by the \code{wilcoxPairedR} function with one group equal
#'           to a vector of \code{mu}.
#'           The author knows of no reference for this technique.
#'           
#'           Currently, the function makes no provisions for \code{NA}
#'           values in the data.  It is recommended that \code{NA}s be removed
#'           beforehand.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_02.html}
#' @concept effect size ordinal nominal
#' @return A single statistic, r
#'         
#' @examples
#' data(Pooh)
#' Data = Pooh[Pooh$Time==2,]
#' wilcox.test(Data$Likert, mu=3, exact=FALSE)
#' wilcoxonOneSampleR(x = Data$Likert, mu=3)
#' 
#' @importFrom coin wilcoxsign_test
#' 
#' @export
 
wilcoxonOneSampleR = function (x, mu=NULL, digits=3, ... ){

  n=length(x)
  
  MU = rep(mu, n)
    
  WT = suppressWarnings(wilcoxsign_test(x ~ MU, ...))
  
  Z = as.numeric(statistic(WT, type="standardized"))
  r = abs(Z)/sqrt(n)
  names(r) = "r"
  R = signif(r, digits=digits)
  
  return(R)
}
