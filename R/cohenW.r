#' @title Cohen's w (omega)
#'
#' @description Calculates Cohen's w for a table of nominal variables.
#' 
#' @param x Either a two-way table or a two-way matrix.
#'          Can also be a vector of observations for one dimension
#'          of a two-way table. 
#' @param y If \code{x} is a vector, \code{y} is the vector of observations for
#'          the second dimension of a two-way table.
#' @param p If \code{x} is a vector of observed counts, \code{p} can be given as
#'          a vector of expected probabilties,
#'          as in a chi-square goodness of fit test.
#' @param digits The number of significant digits in the output.
#' @param ...    Additional arguments passed to \code{chisq.test}.             
#' 
#' @details  Cohen's w is used as a measure of association
#'           between two nominal variables, or as an effect size
#'           for a chi-square test of association.  For a 2 x 2 table,
#'           the absolute value of the phi statistic is the same as
#'           Cohen's w.  
#'           The value of Cohen's w is not bound by 1 on the upper end.
#'           Here, the value is always positive.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/H_10.html}
#' @seealso \code{\link{cramerV}} \code{\link{cramerVFit}}
#' @concept correlation phi cohen w omega
#' @return A single statistic, Cohen's w.
#'         
#' @examples
#' ### Example with table
#' data(Anderson)
#' fisher.test(Anderson)
#' cohenW(Anderson)
#' 
#' ### Example for goodness-of-fit
#' ### Bird foraging example, Handbook of Biological Statistics
#' observed = c(70,   79,   3,    4)
#' expected = c(0.54, 0.40, 0.05, 0.01)
#' chisq.test(observed, p = expected)
#' cohenW(observed, p = expected)
#'
#' ### Example with two vectors
#' Species = c(rep("Species1", 16), rep("Species2", 16))
#' Color   = c(rep(c("blue", "blue", "blue", "green"),4),
#'             rep(c("green", "green", "green", "blue"),4))
#' fisher.test(Species, Color)
#' cohenW(Species, Color)
#' 
#' @importFrom stats chisq.test
#' 
#' @export

cohenW = function(x, y=NULL, p=NULL, digits=4, ...) {
  CW=NULL
  if(is.factor(x)){x=as.vector(x)}
  if(is.factor(y)){x=as.vector(y)}
  if(is.vector(x) & is.vector(y)){
  Chi.sq = suppressWarnings(chisq.test(x, y, correct=FALSE, ...))
  }
  
  if(is.vector(x) & !is.null(p)){
  Chi.sq = suppressWarnings(chisq.test(x=x, p=p, correct=FALSE, ...))
  }
  
 if(is.matrix(x)){x=as.table(x)}
  
 if(is.table(x)){
  Chi.sq = suppressWarnings(chisq.test(x, correct=FALSE, ...))
  }
  
  Sum      = sum(Chi.sq$observed)
  Expected = Chi.sq$expected/Sum
  Observed = Chi.sq$observed/Sum

  CW       = sqrt(sum((Expected-Observed)^2/Expected))
  
  CW = signif(as.numeric(CW), digits=digits)
  names(CW) = "Cohen w"
 return(CW)
}
