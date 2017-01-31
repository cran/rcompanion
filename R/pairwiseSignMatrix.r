#' @title Pairwise sign tests for paired data with matrix output
#'
#' @description Conducts pairwise sign tests across groups
#'              for paired data.

#' @param formula A formula indicating the measurement variable and
#'                the grouping variable. e.g. y ~ group.
#' @param data   The data frame to use.                
#' @param x      The response variable as a vector.
#' @param g      The grouping variable as a vector.
#' @param method The p-value adjustment method to use for multiple tests.
#'               See \code{\link{p.adjust}}.
#' @param ...    Additional arguments passed to
#'               \code{\link{SIGN.test}}.               
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
#' @note    The parsing of the formula is simplistic. 
#'          The first variable on the
#'          left side is used as the measurement variable.  
#'          The first variable on the
#'          right side is used for the grouping variable.         
#'                           
#' @examples
#' data(BobBelcher)
#' friedman.test(Likert ~ Instructor | Rater,
#'               data = BobBelcher)
#' BobBelcher$Instructor = factor( BobBelcher$Instructor, 
#'                   levels = c("Linda Belcher", "Louise Belcher",
#'                              "Tina Belcher", "Bob Belcher",
#'                              "Gene Belcher"))          
#' PT = pairwiseSignMatrix(Likert ~ Instructor,
#'                         data   = BobBelcher,
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