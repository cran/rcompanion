#' @title r effect size for Wilcoxon two-sample paired signed-rank test
#' 
#' @description Calculates r effect size
#'              for a Wilcoxon two-sample paired signed-rank test.
#' 
#' @param x A vector of observations of an ordinal variable.
#' @param g The vector of observations for
#'          the grouping, nominal variable.
#'          Only the first two levels of the nominal variable are used.
#'          The data must be ordered so that the first observation of the
#'          of the first group is paired with the first observation of the
#'          second group.
#' @param digits The number of significant digits in the output.
#' @param ... Additional arguments passed to the \code{wilcoxsign_test} function.             
#' 
#' @details  A Z value is extracted from the \code{wilcoxsign_test} function in the
#'           coin package.  r  is calculated as Z divided by 
#'           square root of the number of observations in one group.
#'
#'           Currently, the function makes no provisions for \code{NA}
#'           values in the data.  It is recommended that \code{NA}s be removed
#'           beforehand.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_06.html}
#' @concept effect size ordinal nominal
#' @return A single statistic, r
#'         
#' @examples
#' data(Pooh)
#' wilcox.test(Likert ~ Time, data=Pooh, paired=TRUE, exact=FALSE)
#' wilcoxonR(x = Pooh$Likert, g = Pooh$Time)
#' 
#' @importFrom coin wilcoxsign_test
#' 
#' @export
 
wilcoxonPairedR = function (x, g=NULL, digits=3, ... ){

  g = factor(g)
  x = x[as.numeric(g)<3]
  g = g[as.numeric(g)<3]
  g = droplevels(g)
    
  WT = suppressWarnings(wilcoxsign_test(x[as.numeric(g)==1] ~ x[as.numeric(g)==2], ...))
  
  Z = as.numeric(statistic(WT, type="standardized"))
  n = length(x[as.numeric(g)==1])
  r = abs(Z)/sqrt(n)
  names(r) = "r"
  R = signif(r, digits=digits)
  
  return(R)
}
