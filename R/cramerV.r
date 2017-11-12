#' @title Cramer's V (phi)
#'
#' @description Calculates Cramer's V for a table of nominal variables.
#' 
#' @param x Either a two-way table or a two-way matrix.
#'          Can also be a vector of observations for one dimension
#'          of a two-way table. 
#' @param y If \code{x} is a vector, \code{y} is the vector of observations for
#'          the second dimension of a two-way table.
#' @param digits The number of significant digits in the output.              
#' @param bias.correct If \code{TRUE}, a bias correction is applied.
#' @param ...    Additional arguments passed to \code{chisq.test}. 
#' 
#' @details  Cramer's V is used as a measure of association
#'           between two nominal variables, or as an effect size
#'           for a chi-square test of association.  For a 2 x 2 table,
#'           the absolute value of the phi statistic is the same as
#'           Cramer's V.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/H_10.html}
#' @concept correlation phi cramer V
#' @seealso \code{\link{cohenW}}
#' @return A single statistic, Cramer's V.
#'         
#' @examples
#' ### Example with table
#' data(Anderson)
#' fisher.test(Anderson)
#' cramerV(Anderson)
#'
#' ### Example with two vectors
#' Species = c(rep("Species1", 16), rep("Species2", 16))
#' Color   = c(rep(c("blue", "blue", "blue", "green"),4),
#'             rep(c("green", "green", "green", "blue"),4))
#' fisher.test(Species, Color)
#' cramerV(Species, Color)
#' 
#' @importFrom stats chisq.test
#' 
#' @export

cramerV = function(x, y=NULL, digits=4, bias.correct=FALSE, ...) {
  CV=NULL
  if(is.factor(x)){x=as.vector(x)}
  if(is.factor(y)){x=as.vector(y)}
  if(is.vector(x) & is.vector(y)){
  N = length(x)
  Chi.sq = suppressWarnings(chisq.test(x, y, correct=FALSE, ...)$statistic)
  Phi =  Chi.sq / N
  R   = length(unique(x))
  C   = length(unique(y))
  CV =  sqrt(Phi / min(R-1, C-1))
  }
  
 if(is.matrix(x)){x=as.table(x)}
  
 if(is.table(x)){
  N = sum(x)
  Chi.sq = suppressWarnings(chisq.test(x, correct=FALSE, ...)$statistic)
  Phi =  Chi.sq / N
  R   = nrow(x)
  C   = ncol(x)
  CV =  sqrt(Phi / min(R-1, C-1))}

 if(bias.correct==TRUE){Phi = max(0, Phi-((R-1)*(C-1)/(N-1)))
                        CC  = C-((C-1)^2/(N-1))
                        RR  = R-((R-1)^2/(N-1))
                        CV  = sqrt(Phi / min(RR-1, CC-1))}  
  
  CV = signif(as.numeric(CV), digits=digits)
  names(CV) = "Cramer V"
 return(CV)
}
