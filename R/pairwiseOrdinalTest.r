#' @title Pairwise two-sample ordinal regression
#'
#' @description Performs pairwise two-sample ordinal regression across groups.
#' 
#' @param x      The response variable as a vector.
#' @param g      The grouping variable as a vector.
#' @param method The p-value adjustment method to use for multiple tests.
#'               See \code{\link{p.adjust}}.
#' @param ...    Additional arguments passed to
#'               \code{\link{clm}}.               
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
#'          The \code{pairwiseOrdinalTest} function 
#'          can be used as a post-hoc method following an omnibus 
#'          ordinal regession whose form is analogous to
#'          a one-way analysis of variance.
#'          
#'          The \code{x} variable must be an ordered factor.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/G_07.html}
#' @seealso \code{\link{pairwiseOrdinalMatrix}}
#' @concept ordinal regression post-hoc one-way
#' @return A dataframe of the groups being compared, the p-values,
#'         and the adjusted p-values. 
#'         
#' @examples
#' data(PoohPiglet)
#' PoohPiglet$Likert.f = factor(PoohPiglet$Likert, ordered = TRUE)
#' PoohPiglet$Speaker = factor(PoohPiglet$Speaker, 
#'                         levels=c("Pooh", "Tigger", "Piglet"))
#' PT = pairwiseOrdinalTest(x      = PoohPiglet$Likert.f,
#'                          g      = PoohPiglet$Speaker,
#'                          method = "fdr")
#' PT
#' cldList(comparison = PT$Comparison,
#'         p.value    = PT$p.adjust,
#'         threshold  = 0.05)
#'  
#' @importFrom stats p.adjust anova
#' @importFrom ordinal clm
#' 
#' @export

pairwiseOrdinalTest = 
  function(x, g, method = "fdr", ...)
  {
  if(!is.factor(g)){g=factor(g)}
  n = length(levels(g))
  N = n*(n-1)/2
  d = data.frame(x = x, g = g)
  Z = data.frame(Comparison=rep("A", N),
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
     z = clm(x ~ g2, data=Dataz, ...)
     y = clm(x ~ 1 , data=Dataz, ...)
     P = signif(anova(z,y)$"Pr(>Chisq)"[2], digits=4)
     P.adjust = NA                       
     Z[k,] =c( paste0(Namea, " - ", Nameb, " = 0"), 
             P, P.adjust)
     }
    } 
  Z$p.adjust = signif(p.adjust(Z$p.value, method = method), digits=4) 
  return(Z)
  }