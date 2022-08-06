#' @title Mangiafico's d
#'
#' @description Calculates Mangiafico's d,
#'              which is the difference in medians
#'              divided by the pooled median absolute deviation,
#'              for several groups in a pairwise manner.
#' 
#' @param formula A formula indicating the response variable and
#'                the independent variable. e.g. y ~ group.
#' @param data   The data frame to use. 
#' @param x If no formula is given, the response variable.
#' @param g If no formula is given, the grouping variable.
#' @param digits The number of significant digits in the output.
#' @param ... Additional arguments passed to the \code{mad()} function. 
#'             
#' @details Mangiafico's d is an appropriate effect size statistic where
#'          Mood's median test, or another test comparing two medians,
#'          might be used.  Note that the response variable is treated
#'          as at least interval.
#'          
#'          When the data in the first group are greater than
#'          in the second group, d is positive.
#'          When the data in the second group are greater than
#'          in the first group, d is negative.
#'          
#'           Be cautious with this interpretation, as R will alphabetize
#'           groups in the formula interface if the grouping variable
#'           is not already a factor.
#'          
#'           Currently, the function makes no provisions for \code{NA}
#'           values in the data.  It is recommended that \code{NA}s be removed
#'           beforehand.
#'                      
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_09.html}
#' @seealso \code{\link{mangiaficoD}}
#' @concept effect size
#' @return A list containing a data frame of pairwise statistics,
#'         and the comparison with the most extreme value
#'         of the statistic.
#'         
#' @note    The parsing of the formula is simplistic. 
#'          The first variable on the
#'          left side is used as the measurement variable.  
#'          The first variable on the
#'          right side is used for the grouping variable.   
#'          
#' @examples
#' data(Catbus)
#' multiMangiaficoD(Steps ~ Teacher, data=Catbus)
#' 
#' @importFrom stats median mad
#' @importFrom boot boot boot.ci
#' 
#' @export
#'
multiMangiaficoD = function(formula=NULL, data=NULL, 
    x=NULL, g=NULL,
    digits = 3, ...)
  {
  if(!is.null(formula)){
    x  = eval(parse(text=paste0("data","$",all.vars(formula[[2]])[1])))
    g  = eval(parse(text=paste0("data","$",all.vars(formula[[3]])[1])))
    }
  if(!is.factor(g)){g=factor(g)}
  l = length(levels(g))
  N = l*(l-1)/2
  d = data.frame(x = x, g = g)
  
  W = data.frame(Comparison   = rep("A", N),
                 Median.1     = rep(NA, N),
                 Median.2     = rep(NA, N),
                 MAD.1        = rep(NA, N),
                 MAD.2        = rep(NA, N),
                 Difference   = rep(NA, N),
                 Pooled.MAD   = rep(NA, N),
                 Mangiafico.d = rep(NA, N),
                 stringsAsFactors=FALSE)
  
  k=0
  for(i in 1:(l-1)){
     for(j in (i+1):l){
     k=k+1
     Namea = as.character(levels(g)[i])
     Nameb = as.character(levels(g)[j])
     Dataz = rbind(subset(d, g==levels(g)[i]), subset(d, g==levels(g)[j]))
     Dataz$g2 = factor(Dataz$g)
     
     m1         = median(Dataz$x[Dataz$g2==levels(g)[i]])
     m2         = median(Dataz$x[Dataz$g2==levels(g)[j]])
     mad1       = mad(Dataz$x[Dataz$g2==levels(g)[i]], ...)
     mad2       = mad(Dataz$x[Dataz$g2==levels(g)[j]], ...)
     difference = signif((m1 - m2), digits=digits)
     pooled     = sqrt((mad1^2 + mad2^2)/2)
     mangia     = (difference / pooled) / pooled
     
     M1     = signif(m1, digits=digits)
     M2     = signif(m2, digits=digits)
     MAD1   = signif(mad1, digits=digits)
     MAD2   = signif(mad2, digits=digits)
     DIFF   = signif(difference, digits=digits)
     POOL   = signif(pooled, digits=digits)
     MANGIA = signif(mangia, digits=digits)
     
     W[k,1] = paste0(Namea, " - ", Nameb)
     
     W[k,2:8] = c(M1, M2, MAD1, MAD2, DIFF, POOL, MANGIA)
     
     }
  }
  
    STAT.m = signif(max(abs(W$Mangiafico.d)), digits=digits)
    COMP   = W$Comparison[which(abs(W$Mangiafico.d)==STAT.m)]
    names(STAT.m) = "Maximum(abs(d))"
    names(COMP)   = rep("Comparison", length(COMP))
  
  Z = list(pairs = W,
           comparison = COMP,
           statistic.m = STAT.m)
    
  return(Z)
  }
