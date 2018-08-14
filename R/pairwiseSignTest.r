#' @title Pairwise sign tests for paired data
#'
#' @description Conducts pairwise sign tests across groups
#'              for paired data.
#'              
#' @param formula A formula indicating the measurement variable and
#'                the grouping variable. e.g. y ~ group.
#' @param data   The data frame to use.
#' @param x      The response variable as a vector.
#' @param g      The grouping variable as a vector.
#' @param method The p-value adjustment method to use for multiple tests.
#'               See \code{stats::p.adjust}.
#' @param ...    Additional arguments passed to
#'               \code{BSDA::SIGN.test}.              
#'             
#' @details The input should include either \code{formula} and \code{data};
#'          or \code{x}, and \code{g}.
#'          
#'          The two sample paired sign test compares medians 
#'          among two groups with paired data.
#'          See \url{http://rcompanion.org/handbook/F_07.html} for
#'          futher discussion of this test.
#' 
#'          The \code{pairwiseSignTest} function
#'          can be used as a post-hoc method following an omnibus
#'          Friedman test.
#'          
#'          The function assumes that the data frame is already ordered by
#'          the blocking variable, so that the first observation of Group 1
#'          is paired with the first observation of Group 2, and so on.                                                                                              
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_10.html}
#' @seealso \code{\link{pairwiseSignMatrix}}
#' @concept median nonparametric post-hoc paired Friedman unreplicated
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
#' data(BobBelcher)
#' BobBelcher = BobBelcher[order(BobBelcher$Instructor, BobBelcher$Rater),]
#' friedman.test(Likert ~ Instructor | Rater,
#'               data = BobBelcher)
#' BobBelcher$Instructor = factor( BobBelcher$Instructor, 
#'                   levels = c("Linda Belcher", "Louise Belcher",
#'                              "Tina Belcher", "Bob Belcher",
#'                              "Gene Belcher"))
#' PT = pairwiseSignTest(Likert ~ Instructor,
#'                       data   = BobBelcher,
#'                       method = "fdr")
#' PT
#' cldList(comparison = PT$Comparison,
#'         p.value    = PT$p.adjust,
#'         threshold  = 0.05)
#'         
#' @importFrom stats p.adjust
#' @importFrom BSDA SIGN.test
#' 
#' @export

pairwiseSignTest = 
  function(formula=NULL, data=NULL,
           x=NULL, g=NULL, 
           method = "fdr", ...)
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
                 S=rep(NA, N),
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
     z = SIGN.test(Datax$x, Datay$x, conf.level=1, ...)
     P = signif(z$p.value, digits=4)
     S = signif(z$statistic, digits=4)
     P.adjust = NA                       
     Z[k,] =c( paste0(Namea, " - ", Nameb, " = 0"), 
             S, P, P.adjust)
     }
    } 
  Z$p.adjust = signif(p.adjust(Z$p.value, method = method), digits=4) 
  return(Z)
  }