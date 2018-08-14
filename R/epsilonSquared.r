#' @title Epsilon-squared
#'
#' @description Calculates epsilon-squared for a table with one ordinal
#'              variable and one nominal variable.
#' 
#' @param x Either a two-way table or a two-way matrix.
#'          Can also be a vector of observations of an ordinal variable.
#' @param g If \code{x} is a vector, \code{g} is the vector of observations for
#'          the grouping, nominal variable.
#' @param group If \code{x} is a table or matrix, \code{group} indicates whether
#'              the \code{"row"} or the \code{"column"} variable is
#'              the nominal, grouping variable. 
#' @param digits The number of significant digits in the output.
#' @param ... Additional arguments passed to the \code{kruskal.test} function.             
#' 
#' @details  Epsilon-squared is used as a measure of association
#'           for the Kruskal-Wallis test or for a two-way
#'           table with one ordinal and one nominal variable.
#'
#'           Currently, the function makes no provisions for \code{NA}
#'           values in the data.  It is recommended that \code{NA}s be removed
#'           beforehand.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/H_11.html}
#' @seealso \code{\link{freemanTheta}}
#' @concept correlation epsilon ordinal nominal
#' @return A single statistic, epsilon-squared
#'         
#' @examples
#' data(Breakfast)
#' library(coin)
#' chisq_test(Breakfast, scores = list("Breakfast" = c(-2, -1, 0, 1, 2)))
#' epsilonSquared(Breakfast)
#' 
#' data(PoohPiglet)
#' kruskal.test(Likert ~ Speaker, data = PoohPiglet)
#' epsilonSquared(x = PoohPiglet$Likert, g = PoohPiglet$Speaker)
#' 
#' ### Same data, as matrix of counts
#' data(PoohPiglet)
#' XT = xtabs( ~ Speaker + Likert , data = PoohPiglet)
#' epsilonSquared(XT)
#' 
#' @importFrom stats kruskal.test
#' 
#' @export
 
epsilonSquared = function (x, g=NULL, group="row", digits=3, ... ){
  
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

  g  = factor(g)
  g  = droplevels(g)
  n  = length(g)
  
  KW = kruskal.test(x, g, ...)

  e2 = KW$statistic / (n-1)
  E2 = signif(e2, digits=digits)
  names(E2) = "epsilon.squared"
  
  return(E2)
}
