#' @title Cramer's V for chi-square goodness-of-fit tests
#'
#' @description Calculates Cramer's V for a vector of counts and expected counts.
#' 
#' @param x A vector of observed counts.
#' @param p A vector of expected or default probabilities.
#' @param digits The number of significant digits in the output.              
#' @param ...    Additional arguments passed to \code{chisq.test}. 
#' 
#' @details  In the case of single vector of counts and expected probabilities,
#'           a modification of Cramer's V can be used to indicate 
#'           the degree of deviation
#'           from the expected probabilities.
#'           
#'           It is not affected by sample size and can be used as an effect
#'           size.
#'           
#'           In the case of equally-distributed expected frequencies,
#'           Cramer's V will be equal to 1 when all counts are in one category,
#'           and it will be equal to 0 when the counts are equally distributed
#'           across categories.
#' 
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/H_03.html}
#' @concept correlation effect size cramer V
#' @seealso \code{\link{cramerV}}
#' @return A single statistic, Cramer's V.
#'         
#' @examples
#' ### Equal probabilities example
#' ### From http://rcompanion.org/handbook/H_03.html
#' nail.color = c("Red", "None", "White", "Green", "Purple", "Blue")
#' observed   = c( 19,    3,      1,       1,       2,        2    )
#' expected   = c( 1/6,   1/6,    1/6,     1/6,     1/6,      1/6  )
#' chisq.test(x = observed, p = expected)
#' cramerVFit(x = observed, p = expected)
#' 
#' ### Unequal probabilities example
#' ### From http://rcompanion.org/handbook/H_03.html
#' race = c("White", "Black", "American Indian", "Asian", "Pacific Islander",
#'           "Two or more races")
#' observed = c(20, 9, 9, 1, 1, 1)
#' expected = c(0.775, 0.132, 0.012, 0.054, 0.002, 0.025)
#' chisq.test(x = observed, p = expected)
#' cramerVFit(x = observed, p = expected)
#' 
#' ### Examples of perfect and zero fits
#' cramerVFit(c(100, 0, 0, 0, 0))
#' cramerVFit(c(10, 10, 10, 10, 10))
#' 
#' @importFrom stats chisq.test
#' 
#' @export

cramerVFit = function(x, p=rep(1/length(x), length(x)), digits=4, ...) {
  CV=NULL
  N = sum(x)
  Chi.sq = suppressWarnings(chisq.test(x=x, p=p, ...)$statistic)
  K   = length(x)
  CV =  sqrt(Chi.sq/N/(K-1))
  
  CV = signif(as.numeric(CV), digits=digits)
  names(CV) = "Cramer V"
 return(CV)
}
