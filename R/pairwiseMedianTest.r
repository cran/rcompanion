#' @title Pairwise Mood's median tests
#'
#' @description Conducts pairwise Mood's median tests across groups.
#' 
#' @param formula A formula indicating the measurement variable and
#'                the grouping variable. e.g. y ~ group.
#' @param data   The data frame to use.
#' @param x      The response variable as a vector.
#' @param g      The grouping variable as a vector.
#' @param exact  If \code{TRUE}, then asks the \code{mood.medtest} function
#'               to conduct an exact test. If \code{NULL}, then uses an 
#'               exact test if the number of values is less than 200.
#'               See \code{\link{mood.medtest}}.
#' @param method The p-value adjustment method to use for multiple tests.
#'               See \code{\link{p.adjust}}.
#' @param ...    Additional arguments passed to
#'               code{\link{mood.medtest}}.               
#'             
#' @details The input should include either \code{formula} and \code{data};
#'          or \code{x}, and \code{g}.
#'          
#'          Mood's median test compares medians among two or more groups.
#'          See \url{http://rcompanion.org/handbook/F_09.html} for
#'          futher discussion of this test.
#' 
#'          The \code{pairwiseMedianTest} function
#'          can be used as a post-hoc method following an omnibus Mood's
#'          median test.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_09.html}
#' @seealso \code{\link{pairwiseMedianMatrix}}
#' @concept moods median nonparametric post-hoc one-way
#' @return A dataframe of the groups being compared, the p-values,
#'         and the adjusted p-values. 
#'         
#' @note    The parsing of the formula is simplistic. 
#'          The first variable on the
#'          left side is used as the measurement variable.  
#'          The first variable on the
#'          right side is used for the grouping variable.   
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
#' @importFrom RVAideMemoire mood.medtest
#' 
#' @export

pairwiseMedianTest = 
  function(formula=NULL, data=NULL, 
    x=NULL, g=NULL, 
    exact = NULL, method = "fdr", ...)
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
     z = mood.medtest(x=Dataz$x, g=Dataz$g2, exact=exact, ...)
     P = signif(z$p.value, digits=4)
     P.adjust = NA                       
     Z[k,] =c( paste0(Namea, " - ", Nameb, " = 0"), 
             P, P.adjust)
     }
    } 
  Z$p.adjust = signif(p.adjust(Z$p.value, method = method), digits=4) 
  return(Z)
  }