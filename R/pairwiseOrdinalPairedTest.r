#' @title Pairwise two-sample ordinal regression for paired or blocked data
#'
#' @description Performs pairwise two-sample ordinal regression across groups
#'              for paired or blocked data.
#' 
#' @param x      The response variable as a vector.
#' @param g      The grouping variable as a vector.
#' @param b      The blocking variable as a vector.
#' @param method The p-value adjustment method to use for multiple tests.
#'               See \code{\link{p.adjust}}.
#' @param ...    Additional arguments passed to
#'               \code{\link{clmm}}.               
#'             
#' @details 
#'          Ordinal regression 
#'          is analogous to general linear regression or generalized linear
#'          regression for cases where 
#'          the dependent variable
#'          is an ordinal variable.
#'          The \code{ordinal} package provides a flexible and powerful
#'          implementation of ordinal regression.
#'          
#'          The \code{pairwiseOrdinalPairedTest} function 
#'          can be used as a post-hoc method following an omnibus 
#'          ordinal regession whose form is analogous to
#'          a one-way analysis of variance with random blocks.
#'          
#'          The blocking variable is treated as a random variable.
#'          
#'          The \code{x} variable must be an ordered factor.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/G_08.html}
#' @seealso \code{\link{pairwiseOrdinalPairedMatrix}}
#' @concept ordinal regression post-hoc one-way paired blocks
#' @return A dataframe of the groups being compared, the p-values,
#'         and the adjusted p-values. 
#'         
#' @examples
#' data(BobBelcher)
#' BobBelcher$Likert.f = factor(BobBelcher$Likert, ordered = TRUE)
#' BobBelcher = BobBelcher[order(factor(BobBelcher$Instructor, 
#'                         levels=c("Linda Belcher", "Louise Belcher",
#'                                  "Tina Belcher", "Bob Belcher",
#'                                  "Gene Belcher"))),] 
#' pairwiseOrdinalPairedTest(x      = BobBelcher$Likert.f,
#'                           g      = BobBelcher$Instructor,
#'                           b      = BobBelcher$Rater,
#'                           threshold="equidistant",
#'                           method = "fdr")
#' 
#' @importFrom stats p.adjust pchisq
#' @importFrom ordinal clmm
#' 
#' @export

pairwiseOrdinalPairedTest = 
  function(x, g, b, method = "fdr", ...)
  {
  n = length(unique(g))
  N = n*(n-1)/2
  d = data.frame(x = x, g = g, b = b)
  Z = data.frame(Comparison=rep("A", N),
                 p.value=rep(NA, N),
                 p.adjust=rep(NA, N),
                 stringsAsFactors=FALSE) 
  k=0               
  for(i in 1:(n-1)){
     for(j in (i+1):n){
       k=k+1
     Namea = as.character(unique(g)[i])
     Nameb = as.character(unique(g)[j])
     Datax = subset(d, g==unique(g)[i])
     Datay = subset(d, g==unique(g)[j])
     Dataz = rbind(Datax, Datay)
     Dataz$g2 = factor(Dataz$g)
     z = clmm(x ~ g2 + (1|b), data=Dataz, ...)
     y = clmm(x ~  1 + (1|b), data=Dataz, ...)
     p = pchisq(2*(z$logLik-y$logLik),1,lower.tail=FALSE)
     P = signif(p, digits=4)
     P.adjust = NA                       
     Z[k,] =c( paste0(Namea, " - ", Nameb, " = 0"), 
             P, P.adjust)
     }
    } 
  Z$p.adjust = signif(p.adjust(Z$p.value, method = method), digits=4) 
  return(Z)
  }