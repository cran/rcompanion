#' @title Pairwise two-sample ordinal regression for blocked data with matrix
#'        output      
#'
#' @description Performs pairwise two-sample ordinal regression across groups
#'              for paired or blocked data.
#' @param formula A formula indicating the measurement variable and
#'                the grouping variable. e.g. y ~ group | block.
#' @param data   The data frame to use.                
#' @param x      The response variable as a vector.
#' @param g      The grouping variable as a vector.
#' @param b      The blocking variable as a vector.
#' @param method The p-value adjustment method to use for multiple tests.
#'               See \code{\link{p.adjust}}.
#' @param ...    Additional arguments passed to
#'               \code{\link{clmm}}.               
#'             
#' @details The input should include either \code{formula} and \code{data};
#'          or \code{x}, \code{g}, and \code{b}.
#'          
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
#'          The matrix output can be converted to a compact letter display.
#'          
#'          The blocking variable is treated as a random variable.
#'          
#'          The \code{x} variable must be an ordered factor.                                  
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/G_08.html}
#' @seealso \code{\link{pairwiseOrdinalPairedTest}}
#' @concept ordinal regression post-hoc one-way paired blocks cld
#' @return A list consisting of:
#'         A matrix of p-values;
#'         The p-value adjustment method;
#'         A matrix of adjusted p-values. 
#'
#' @note    The parsing of the formula is simplistic. 
#'          The first variable on the
#'          left side is used as the measurement variable.  
#'          The first variable on the
#'          right side is used for the grouping variable.
#'          The second variable on the
#'          right side is used for the blocking variable
#'
#' @examples
#' data(BobBelcher)
#' BobBelcher$Likert.f = factor(BobBelcher$Likert, ordered = TRUE)
#' BobBelcher$Instructor = factor( BobBelcher$Instructor, 
#'                   levels = c("Linda Belcher", "Louise Belcher",
#'                              "Tina Belcher", "Bob Belcher",
#'                              "Gene Belcher")) 
#' PT = pairwiseOrdinalPairedMatrix(Likert.f ~ Instructor | Rater,
#'                                  data      = BobBelcher,
#'                                  threshold ="equidistant",
#'                                  method    = "fdr")$Adjusted
#' PT
#' library(multcompView)
#' multcompLetters(PT,
#'                 compare="<",
#'                 threshold=0.05,
#'                 Letters=letters)
#'                  
#' @importFrom stats p.adjust pchisq
#' @importFrom ordinal clmm
#' 
#' @export

pairwiseOrdinalPairedMatrix = 
  function(formula=NULL, data=NULL,
           x=NULL, g=NULL, b=NULL, 
           method = "fdr", ...)
  {
  if(!is.null(formula)){
    x  = eval(parse(text=paste0("data","$",all.vars(formula[[2]])[1])))
    g  = eval(parse(text=paste0("data","$",all.vars(formula[[3]])[1])))
    b  = eval(parse(text=paste0("data","$",all.vars(formula[[3]])[2])))
    }
  n = length(levels(g))
  N = n*n
  d = data.frame(x = x, g = g, b = b)
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
     z = clmm(x ~ g2 + (1|b), data=Dataz, ...)
     y = clmm(x ~  1 + (1|b), data=Dataz, ...)
     p = pchisq(2*(z$logLik-y$logLik),1,lower.tail=FALSE)           
   Y[i,j] = signif(p, digits=4)
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