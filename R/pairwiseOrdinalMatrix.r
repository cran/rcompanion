#' @title Pairwise two-sample ordinal regression with matrix output
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
#'          The \code{pairwiseOrdinalMatrix} function 
#'          can be used as a post-hoc method following an omnibus 
#'          ordinal regession whose form is analogous to
#'          a one-way analysis of variance.
#'          The matrix output can be converted to a compact letter display.
#'          
#'          The \code{x} variable must be an ordered factor.                                   
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/G_07.html}
#' @seealso \code{\link{pairwiseOrdinalTest}}
#' @concept ordinal regression post-hoc one-way cld
#' @return A list consisting of:
#'         A matrix of p-values;
#'         the p-value adjustment method;
#'         a matrix of adjusted p-values. 
#'         
#' @examples
#' data(PoohPiglet)
#' PoohPiglet$Likert.f = factor(PoohPiglet$Likert, ordered = TRUE)
#' PoohPiglet = PoohPiglet[order(factor(PoohPiglet$Speaker, 
#'                         levels=c("Pooh", "Tigger", "Piglet"))),]               
#' PT = pairwiseOrdinalMatrix(x      = PoohPiglet$Likert.f,
#'                            g      = PoohPiglet$Speaker,
#'                            method = "fdr")$Adjusted
#' PT                          
#' library(multcompView)
#' multcompLetters(PT,
#'                 compare="<",
#'                 threshold=0.05,
#'                 Letters=letters)
#'                  
#' @importFrom stats p.adjust anova
#' @importFrom ordinal clm
#' 
#' @export  

pairwiseOrdinalMatrix = 
  function(x, g, method = "fdr", ...)
  {
  n = length(unique(g))
  N = n*n
  d = data.frame(x = x, g = g)
  Y = matrix(rep(NA_real_, N),ncol=n)
  rownames(Y)=unique(g)
  colnames(Y)=unique(g)
  Z = matrix(rep(NA_real_, N),ncol=n)
  rownames(Z)=unique(g)
  colnames(Z)=unique(g)
  k=0
  for(i in 1:(n-1)){
     for(j in (i+1):n){
     k=k+1
     Datax = subset(d, g==unique(g)[i])
     Datay = subset(d, g==unique(g)[j])
     Dataz = rbind(Datax, Datay)
     Dataz$g2 = factor(Dataz$g)
     z = clm(x ~ g2, data=Dataz, ...)
     y = clm(x ~ 1 , data=Dataz, ...)
   Y[i,j] = signif(anova(z,y)$"Pr(>Chisq)"[2], digits=4)
   } 
   }
Z[upper.tri(Z)] = 
      signif(p.adjust(Y[upper.tri(Y)], method=method), digits=4)
Z = t(Z)
Z[upper.tri(Z)] = 
      signif(p.adjust(Y[upper.tri(Y)], method=method), digits=4)
diag(Z) = signif(1.00, digits = 4)
W = method
V = list(Y, W, Z)
names(V) = c("Unadjusted",
             "Method",
             "Adjusted")
return(V)   
} 