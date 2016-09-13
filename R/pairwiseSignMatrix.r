#' @title Pairwise sign tests for paired data with matrix output
#'
#' @description Conducts pairwise sign tests across groups
#'              for paired data.
#' 
#' @param x      The response variable as a vector.
#' @param g      The grouping variable as a vector.
#' @param method The p-value adjustment method to use for multiple tests.
#'               See \code{\link{p.adjust}}.
#' @param ...    Additional arguments passed to
#'               \code{\link{SIGN.test}}.               
#'             
#' @details The two sample paired sign test compares medians 
#'          among two groups with paired data.
#'          See \url{http://rcompanion.org/handbook/F_07.html} for
#'          futher discussion of this test.
#' 
#'          The \code{pairwiseSignTest} function
#'          can be used as a post-hoc method following an omnibus
#'          Friedman test.
#'          The matrix output can be converted to a compact letter display.
#'          
#'          The function assumes that the data frame is already ordered by
#'          the blocking variable, so that the first observation of Group 1
#'          is paired with the first observation of Group 2, and so on.                                                                                              
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_10.html}
#' @seealso \code{\link{pairwiseSignTest}}
#' @concept median nonparametric post-hoc paired Friedman unreplicated cld
#' @return A list consisting of:
#'         A matrix of p-values;
#'         the p-value adjustment method;
#'         a matrix of adjusted p-values. 
#'         
#' @examples
#' data(BobBelcher)
#' friedman.test(Likert ~ Instructor | Rater,
#'               data = BobBelcher)
#' BobBelcher = BobBelcher[order(factor(BobBelcher$Instructor, 
#'                         levels=c("Linda Belcher", "Louise Belcher",
#'                                  "Tina Belcher", "Bob Belcher",
#'                                  "Gene Belcher"))),]             
#' PT = pairwiseSignMatrix(x      = BobBelcher$Likert,
#'                         g      = BobBelcher$Instructor,
#'                         method = "fdr")$Adjusted
#' PT
#' library(multcompView)
#' multcompLetters(PT,
#'                 compare="<",
#'                 threshold=0.05,
#'                 Letters=letters)
#'                  
#' @importFrom stats p.adjust
#' @importFrom BSDA SIGN.test
#' 
#' @export

  pairwiseSignMatrix = 
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
     z = SIGN.test(Datax$x, Datay$x, conf.level=1, ...)             
   Y[i,j] = signif(z$p.value, digits = 4)
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