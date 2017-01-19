#' @title Pairwise Mood's median tests with matrix output
#'
#' @description Conducts pairwise Mood's median tests across groups.
#' 
#' @param x      The response variable as a vector.
#' @param g      The grouping variable as a vector.
#' @param exact  If \code{TRUE}, then asks the \code{mood.medtest} function
#'               to conduct an exact test. If \code{NULL}, then uses an 
#'               exact test if the number of values is less than 200.
#'               See \code{\link{mood.medtest}}.
#' @param method The p-value adjustment method to use for multiple tests.
#'               See \code{\link{p.adjust}}.
#' @param ...    Additional arguments passed to
#'               \code{\link{mood.medtest}}.               
#'             
#' @details Mood's median test compares medians among two or more groups.
#'          See \url{http://rcompanion.org/handbook/F_09.html} for
#'          futher discussion of this test.
#' 
#'          The \code{pairwiseMedianMatrix} function
#'          can be used as a post-hoc method following an omnibus Mood's
#'          median test.
#'          The matrix output can be converted to a compact letter display.                                                                                    
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_09.html}
#' @seealso \code{\link{pairwiseMedianTest}} 
#' @concept moods median nonparametric post-hoc one-way cld
#' @return A list consisting of:
#'         a matrix of p-values;
#'         the p-value adjustment method;
#'         a matrix of adjusted p-values.
#'         
#' @examples
#' data(PoohPiglet)
#' PoohPiglet$Speaker = factor(PoohPiglet$Speaker, 
#'                             levels=c("Pooh", "Tigger", "Piglet"))               
#' PT = pairwiseMedianMatrix(x      = PoohPiglet$Likert,
#'                           g      = PoohPiglet$Speaker,
#'                           exact  = NULL,
#'                           method = "fdr")$Adjusted
#' PT                           
#' library(multcompView)
#' multcompLetters(PT,
#'                 compare="<",
#'                 threshold=0.05,
#'                 Letters=letters)                    
#' 
#' @importFrom stats p.adjust
#' @importFrom RVAideMemoire mood.medtest
#' 
#' @export

pairwiseMedianMatrix = 
  function(x, g, exact=NULL, method = "fdr", ...)
  {
  if(!is.factor(g)){g=factor(g)}
  n = length(levels(g))
  N = n*n
  d = data.frame(x = x, g = g)
  Y = matrix(rep(NA_real_, N),ncol=n)
  rownames(Y)=levels(g)
  colnames(Y)=levels(g)
  Z = matrix(rep(NA_real_, N),ncol=n)
  rownames(Z)=levels(g)
  colnames(Z)=levels(g)
  k=0
  for(i in 1:(n-1)){
     for(j in (i+1):n){
     k=k+1
     Datax = subset(d, g==levels(g)[i])
     Datay = subset(d, g==levels(g)[j])
     Dataz = rbind(Datax, Datay)
     Dataz$g2 = factor(Dataz$g)
     z = mood.medtest(x=Dataz$x, g=Dataz$g2, exact=exact, ...)             
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