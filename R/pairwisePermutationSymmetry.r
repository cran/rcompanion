#' @title Pairwise two-sample permutation symmetry tests
#'
#' @description Conducts pairwise two-sample permutation tests 
#'              of symmetry across groups.
#' 
#' @param x      The response variable as a vector.
#' @param g      The grouping variable as a vector.
#' @param b      The blocking variable as a vector.
#' @param method The p-value adjustment method to use for multiple tests.
#'               See \code{\link{p.adjust}}.
#' @param ...    Additional arguments passed to
#'               \code{\link{symmetry_test}}.               
#'             
#' @details Permutation tests are non-parametric tests 
#'          that do not assume normally-distributed errors.
#'          See \url{http://rcompanion.org/rcompanion/d_06a.html} for
#'          futher discussion of this test.
#' 
#'          The \code{pairwisePermutationSymmetry} function
#'          can be used as a post-hoc method following an omnibus 
#'          permutation test
#'          analogous to a paired one-way analysis of variance.
#'                                                                                              
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/K_03.html}
#' @seealso \code{\link{pairwisePermutationSymmetryMatrix}}
#' @concept permutation nonparametric post-hoc one-way symmetry
#' @return A dataframe of the groups being compared, the p-values,
#'         and the adjusted p-values. 
#'         
#' @examples
#' data(BobBelcher)
#' BobBelcher$Instructor = factor(BobBelcher$Instructor, 
#'                                levels=c("Linda Belcher", "Louise Belcher",
#'                                  "Tina Belcher", "Bob Belcher",
#'                                  "Gene Belcher"))
#' BobBelcher$Likert.f = factor(BobBelcher$Likert, ordered=TRUE)
#' PT = pairwisePermutationSymmetry(x      = BobBelcher$Likert.f,
#'                                  g      = BobBelcher$Instructor,
#'                                  b      = BobBelcher$Rater,
#'                                  method = "fdr")
#' PT
#' cldList(comparison = PT$Comparison,
#'         p.value    = PT$p.adjust,
#'         threshold  = 0.05)
#' 
#' @importFrom stats p.adjust
#' @importFrom coin symmetry_test statistic
#' 
#' @export

pairwisePermutationSymmetry = 
  function(x, g, b, method = "fdr", ...)
  {
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