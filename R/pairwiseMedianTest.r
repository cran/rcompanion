#' @title Pairwise Mood's median tests
#'
#' @description Conducts pairwise Mood's median tests across groups.
#' 
#' @param formula A formula indicating the measurement variable and
#'                the grouping variable. e.g. y ~ group.
#' @param data   The data frame to use.
#' @param x      The response variable as a vector.
#' @param g      The grouping variable as a vector.
#' @param method The p-value adjustment method to use for multiple tests.
#'               See \code{stats::p.adjust}.
#' @param digits The number of significant digits to round output.
#' @param ...    Additional arguments passed to
#'               \code{coin::median_test}.               
#'             
#' @details The input should include either \code{formula} and \code{data};
#'          or \code{x}, and \code{g}.
#'          
#'          Mood's median test compares medians among two or more groups.
#'          See \url{http://rcompanion.org/handbook/F_09.html} for
#'          further discussion of this test.
#' 
#'          The \code{pairwiseMedianTest} function
#'          can be used as a post-hoc method following an omnibus Mood's
#'          median test.  It passes the data for pairwise groups to
#'          \code{coin::median_test}.
#'          
#'          The output can be converted to a compact letter display,
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
#' @references \url{http://rcompanion.org/handbook/F_09.html}
#' 
#' @seealso \code{\link{pairwiseMedianMatrix}}
#' 
#' @concept post-hoc
#' @concept Mood's median test
#' 
#' @return A dataframe of the groups being compared, the p-values,
#'         and the adjusted p-values. 
#'
#' @examples
#' data(PoohPiglet)
#' PoohPiglet$Speaker = factor(PoohPiglet$Speaker,
#'                      levels = c("Pooh", "Tigger", "Piglet"))
#' PT = pairwiseMedianTest(Likert ~ Speaker,
#'                         data   = PoohPiglet,
#'                         exact  = NULL,
#'                         method = "fdr")
#' PT                         
#' cldList(comparison = PT$Comparison,
#'         p.value    = PT$p.adjust,
#'         threshold  = 0.05)                         
#' 
#' @importFrom stats p.adjust
#' @importFrom coin median_test
#' 
#' @export

pairwiseMedianTest = 
  function(formula=NULL, data=NULL, 
    x=NULL, g=NULL, 
    digits = 4, method = "fdr", ...)
  {
  if(!is.null(formula)){
    x  = eval(parse(text=paste0("data","$",all.vars(formula[[2]])[1])))
    g  = eval(parse(text=paste0("data","$",all.vars(formula[[3]])[1])))
    }
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
     z = median_test(x ~ g2, data=Dataz, ...)
     P = signif(pvalue(z)[1], digits=digits)
     P.adjust = NA                       
     Z[k,] =c( paste0(Namea, " - ", Nameb, " = 0"), 
             P, P.adjust)
     }
    } 
  Z$p.adjust = signif(p.adjust(Z$p.value, method = method), digits=4) 
  return(Z)
  }