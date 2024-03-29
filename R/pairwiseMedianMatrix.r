#' @title Pairwise Mood's median tests with matrix output
#'
#' @description Conducts pairwise Mood's median tests across groups.
#' 
#' @param formula A formula indicating the measurement variable and
#'                the grouping variable. e.g. y ~ group.
#' @param data   The data frame to use.
#' @param x      The response variable as a vector.
#' @param g      The grouping variable as a vector.
#' @param digits The number of significant digits to round output.
#' @param method The p-value adjustment method to use for multiple tests.
#'               See \code{stats::p.adjust}.
#' @param ...    Additional arguments passed to
#'               \code{coin::median_test}.               
#'             
#' @details The input should include either \code{formula} and \code{data};
#'          or \code{x}, and \code{g}.
#'          
#'          Mood's median test compares medians among two or more groups.
#'          See \url{https://rcompanion.org/handbook/F_09.html} for
#'          futher discussion of this test.
#' 
#'          The \code{pairwiseMedianMatrix} function
#'          can be used as a post-hoc method following an omnibus Mood's
#'          median test. It passes the data for pairwise groups to
#'          \code{coin::median_test}.
#'          
#'          The matrix output can be converted to a compact letter display,
#'          as in the example.
#'          
#' @note    The parsing of the formula is simplistic. 
#'          The first variable on the
#'          left side is used as the measurement variable.  
#'          The first variable on the
#'          right side is used for the grouping variable.                                                                                 
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{https://rcompanion.org/handbook/F_09.html}
#' 
#' @seealso \code{\link{pairwiseMedianTest}} 
#' 
#' @concept post-hoc
#' @concept Mood's median test
#' 
#' @return A list consisting of:
#'         a matrix of p-values;
#'         the p-value adjustment method;
#'         a matrix of adjusted p-values.
#'          
#' @examples
#' data(PoohPiglet)
#' PoohPiglet$Speaker = factor(PoohPiglet$Speaker,
#'                           levels = c("Pooh", "Tigger", "Piglet"))
#' PT = pairwiseMedianMatrix(Likert ~ Speaker,
#'                           data   = PoohPiglet,
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
#' @importFrom coin median_test
#' 
#' @export

pairwiseMedianMatrix = 
  function(formula=NULL, data=NULL, 
           x=NULL, g=NULL, digits=4, method = "fdr", ...)
  {
  if(!is.null(formula)){
    x  = eval(parse(text=paste0("data","$",all.vars(formula[[2]])[1])))
    g  = eval(parse(text=paste0("data","$",all.vars(formula[[3]])[1])))
    }
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
     z = median_test(x ~ g2, data=Dataz, ...)             
   Y[i,j] = signif(pvalue(z), digits = digits)
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