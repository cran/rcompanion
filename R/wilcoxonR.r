#' @title r effect size for Wilcoxon two-sample rank-sum test
#' 
#' @description Calculates r effect size
#'              for Mann-Whitney, two-sample rank-sum test,
#'              or a table with an ordinal variable and a
#'              nominal variable with two levels.
#' 
#' @param x Either a two-way table or a two-way matrix.
#'          Can also be a vector of observations of an ordinal variable.
#' @param g If \code{x} is a vector, \code{g} is the vector of observations for
#'          the grouping, nominal variable.
#'          Only the first two levels of the nominal variable are used.
#' @param group If \code{x} is a table or matrix, \code{group} indicates whether
#'              the \code{"row"} or the \code{"column"} variable is
#'              the nominal, grouping variable. 
#' @param digits The number of significant digits in the output.
#' @param ... Additional arguments passed to the \code{wilcox_test} function.             
#' 
#' @details  A Z value is extracted from the \code{wilcox_test} function in the
#'           coin package.  r  is calculated as Z divided by 
#'           square root of the total observations.
#'  
#'           Currently, the function makes no provisions for \code{NA}
#'           values in the data.  It is recommended that \code{NA}s be removed
#'           beforehand.
#'                      
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_04.html}
#' @seealso \code{\link{freemanTheta}}
#' @concept correlation effect size ordinal nominal
#' @return A single statistic, r
#'         
#' @examples
#' data(Breakfast)
#' Table = Breakfast[1:2,]
#' library(coin)
#' chisq_test(Table, scores = list("Breakfast" = c(-2, -1, 0, 1, 2)))
#' wilcoxonR(Table)
#' 
#' data(PoohPiglet)
#' Data = PoohPiglet[PoohPiglet$Speaker %in% c("Pooh", "Piglet"),]
#' wilcox.test(Likert ~ Speaker, data = Data)
#' wilcoxonR(x = Data$Likert, g = Data$Speaker)
#' 
#' @importFrom coin wilcox_test
#' 
#' @export
 
wilcoxonR = function (x, g=NULL, group="row", digits=3, ... ){
  
  if(is.matrix(x)){x=as.table(x)}
  
  if(is.table(x)){
    Counts = as.data.frame(x)
    Long = Counts[rep(row.names(Counts), Counts$Freq), c(1, 2)]
    rownames(Long) = seq(1:nrow(Long))
    if(group=="row"){
       g=factor(Long[,1])
       x=as.numeric(Long[,2])}
    if(group=="column"){
       g=factor(Long[,2])
       x=as.numeric(Long[,1])}
  }

  g = factor(g)
  x = x[as.numeric(g)<3]
  g = g[as.numeric(g)<3]
  g = droplevels(g)
    
  WT = suppressWarnings(wilcox_test(x ~ g, ...))
  
  Z = as.numeric(statistic(WT, type="standardized"))
  N = length(g)
  r = abs(Z)/sqrt(N)
  names(r) = "r"
  R = signif(r, digits=digits)
  
  return(R)
}
