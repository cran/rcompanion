#' @title Pairwise two-sample symmetry tests
#'
#' @description Conducts pairwise two-sample symmetry tests across groups.
#' 
#' @param formula A formula indicating the measurement variable and
#'                the grouping variable. e.g. y ~ group | block.
#' @param data   The data frame to use.                 
#' @param x      The response variable as a vector.
#' @param g      The grouping variable as a vector.
#' @param b      The blocking variable as a vector.
#' @param method The p-value adjustment method to use for multiple tests.
#'               See \code{stats::p.adjust}.
#' @param ...    Additional arguments passed to
#'               \code{coin::symmetry_test}.               
#'             
#' @details The input should include either \code{formula} and \code{data};
#'          or \code{x}, \code{g}, and \code{b}.
#' 
#'          This function is a wrapper for \code{coin::symmetry_test},
#'          passing pairwise groups to the function. It's critical to read
#'          and understand the documentation for this function to understand
#'          its use and options.
#'
#' @note    The parsing of the formula is simplistic. 
#'          The first variable on the
#'          left side is used as the measurement variable.  
#'          The first variable on the
#'          right side is used for the grouping variable.
#'          The second variable on the
#'          right side is used for the blocking variable.
#'                                                                                              
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{http://rcompanion.org/handbook/K_03.html}
#' 
#' @seealso \code{\link{pairwisePermutationSymmetryMatrix}}
#' 
#' @concept post-hoc 
#' @concept permutation test
#' 
#' @return A dataframe of the groups being compared, the p-values,
#'         and the adjusted p-values. 
#'
#' @examples
#' data(BobBelcher)
#' 
#' BobBelcher$Instructor = factor( BobBelcher$Instructor, 
#'                                 levels = c("Linda Belcher", "Louise Belcher",
#'                                            "Tina Belcher", "Bob Belcher",
#'                                            "Gene Belcher"))
#'                                            
#' library(coin)
#' 
#' symmetry_test(Likert ~ Instructor | Rater, data= BobBelcher,
#'               ytrafo   = rank_trafo,
#'               teststat = "quadratic")
#' 
#' PT = pairwisePermutationSymmetry(Likert ~ Instructor | Rater,
#'                                  data     = BobBelcher,
#'                                  ytrafo   = rank_trafo,
#'                                  teststat = "quadratic",
#'                                  method   = "fdr")
#' PT
#' 
#' cldList(comparison = PT$Comparison,
#'         p.value    = PT$p.adjust,
#'        threshold  = 0.05)
#' 
#' @importFrom stats p.adjust
#' @importFrom coin symmetry_test statistic
#' 
#' @export

pairwisePermutationSymmetry = 
  function(formula=NULL, data=NULL,
           x=NULL, g=NULL, b=NULL, 
           method = "fdr", ...)
  {
  if(!is.null(formula)){
    x  = eval(parse(text=paste0("data","$",all.vars(formula[[2]])[1])))
    g  = eval(parse(text=paste0("data","$",all.vars(formula[[3]])[1])))
    b  = eval(parse(text=paste0("data","$",all.vars(formula[[3]])[2])))
    }
  if(!is.factor(g)){g=factor(g)}
  if(!is.factor(b)){g=factor(b)}
  n = length(levels(g))
  N = n*(n-1)/2
  d = data.frame(x = x, g = g, b = b)
  Z = data.frame(Comparison=rep("A", N),
                 Stat=rep(NA, N),
                 p.value=rep(NA, N),
                 p.adjust=rep(NA, N),
                 stringsAsFactors=FALSE)
                 
  k=0               
  for(i in 1:(n-1)){
     for(j in (i+1):n){
       k=k+1
     Namea = as.character(levels(g)[i])
     Nameb = as.character(levels(g)[j])
     Datax = subset(d, g==levels(g)[i])
     Datay = subset(d, g==levels(g)[j])
     Dataz = rbind(Datax, Datay)
     Dataz$g2 = factor(Dataz$g)
     z = symmetry_test(x ~ g2|b, data=Dataz, ...)
     P = signif(pvalue(z), digits=4)
     S = signif(statistic (z), digits=4)
     P.adjust = NA                       
     Z[k,] =c( paste0(Namea, " - ", Nameb, " = 0"), 
             S, P, P.adjust)
     }
    } 
  Z$p.adjust = signif(p.adjust(Z$p.value, method = method), digits=4) 
  return(Z)
  }