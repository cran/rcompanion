#' @title Pairwise Vargha and Delaney's A and Cliff's delta
#'
#' @description Calculates Vargha and Delaney's A (VDA),
#'              Cliff's delta (CD), and r
#'              for several groups in a pairwise manner.
#' 
#' @param formula A formula indicating the response variable and
#'                the independent variable. e.g. y ~ group.
#' @param data   The data frame to use. 
#' @param x If no formula is given, the response variable.
#' @param g If no formula is given, the grouping variable.
#' @param statistic One of \code{"VDA"}, \code{"CD"}, or \code{"r"}.
#'                  This determines which statistic will be
#'                  evaluated to determine the comparison with the
#'                  most divergent groups.
#' @param digits The number of significant digits in the output.
#' @param ... Additional arguments passed to the \code{wilcox.test} function. 
#'             
#' @details VDA and CD are effect size statistic appropriate
#'          in cases where a Wilcoxon-Mann-Whitney test might be used.
#'          Here, the pairwise approach would be used in cases where
#'          a Kruskal-Wallis test might be used.
#'          VDA ranges from 0 to 1, 
#'          with 0.5 indicating stochastic equality,
#'          and 1 indicating that the first group dominates the second.
#'          CD ranges from -1 to 1, with 0 indicating stochastic equality,
#'          and 1 indicating that the first group dominates the second.
#'          r ranges from approximately, -0.86 to 0.86, 
#'          depending on sample size,
#'          with 0 indicating no effect,
#'          and a positive result indicating
#'          that values in the first group are greater than in the second.
#'          
#'          In the function output,
#'          \code{VDA.m} is the greater of VDA or 1-VDA.
#'          \code{CD.m} is the absolute value of CD.
#'          \code{r.m} is the absolute value of r.
#'          
#'          The function calculates VDA and Cliff's delta from the "W" 
#'          U statistic from the
#'          \code{wilcox.test} function.
#'          Specifically, \code{VDA = U/(n1*n2); CD = (VDA-0.5)*2}.
#'          For r, the Z value is extracted 
#'          from the \code{wilcox_test} function in the
#'           coin package.  r  is calculated as Z divided by 
#'           square root of the total observations.
#'            
#'          The input should include either \code{formula} and \code{data};
#'          or \code{var}, and \code{group}.
#'          
#'           Currently, the function makes no provisions for \code{NA}
#'           values in the data.  It is recommended that \code{NA}s be removed
#'           beforehand.
#'           
#'           When the data in the first group are greater than
#'           in the second group, 
#'           VDA is > 0.5, CD is positive, and r is positive.
#'           When the data in the second group are greater than
#'           in the first group, 
#'           VDA is < 0.5, CD is negative, and r is negative.
#'           Be cautious with this interpretation, as R will alphabetize
#'           groups in the formula interface if the grouping variable
#'           is not already a factor.
#'                      
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_08.html}
#' @seealso \code{\link{cliffDelta}}, \code{\link{cliffDelta}}
#' @concept effect size
#' @return A list containing a data frame of pairwise statistics,
#'         and the comparison with the most extreme value
#'         of the chosen statistic.
#'         
#' @note    The parsing of the formula is simplistic. 
#'          The first variable on the
#'          left side is used as the measurement variable.  
#'          The first variable on the
#'          right side is used for the grouping variable.   
#'          
#' @examples
#' data(PoohPiglet)
#' multiVDA(Likert ~ Speaker, data=PoohPiglet)
#' 
#' @importFrom stats wilcox.test
#' @importFrom boot boot boot.ci
#' 
#' @export
#'
multiVDA = function(formula=NULL, data=NULL, 
    x=NULL, g=NULL, statistic="VDA",
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
  
  W = data.frame(Comparison        = rep("A", N),
                 VDA               = rep(NA, N),
                 CD                = rep(NA, N),
                 r                 = rep(NA, N),
                 VDA.m             = rep(NA, N),
                 CD.m              = rep(NA, N),
                 r.m               = rep(NA, N),
                 stringsAsFactors=FALSE)
  
  k=0
  for(i in 1:(l-1)){
     for(j in (i+1):l){
     k=k+1
     Namea = as.character(levels(g)[i])
     Nameb = as.character(levels(g)[j])
     Dataz = rbind(subset(d, g==levels(g)[i]), subset(d, g==levels(g)[j]))
     Dataz$g2 = factor(Dataz$g)
     
     Test = suppressWarnings(wilcox_test(x ~ g2, data=Dataz))
     Z = as.numeric(statistic(Test, type="standardized"))
     R = signif(Z/sqrt(length(Dataz$g2)), digits=digits)
   
     a = Dataz$x[Dataz$g2==levels(g)[i]]
     b = Dataz$x[Dataz$g2==levels(g)[j]]
     m = as.numeric(length(a))
     n = as.numeric(length(b))
     WT = suppressWarnings(wilcox.test(a, b), ...)$statistic
     VDA = signif(WT/(m*n), digits=digits)

     CD = signif((VDA * 2) - 1, digits=digits)
                   
     R.m = abs(R)
     VDA.m = max(c(VDA, 1-VDA))
     CD.m = abs(CD)
     
     W[k,1] = paste0(Namea, " - ", Nameb)
     
     W[k,2:7] = c(VDA, CD, R, VDA.m, CD.m, R.m)
     
     }
  }
  
  if(statistic == "VDA"){
    STAT   = W$VDA[which(W$VDA.m==max(W$VDA.m))]
    STAT.m = max(W$VDA.m)
    COMP   = W$Comparison[which(W$VDA.m==max(W$VDA.m))]
    names(STAT)   = "VDA"
    names(STAT.m) = "VDA.m"
    names(COMP)   = "Comparison"
      }
  
  if(statistic == "CD"){
    STAT   = W$CD[which(W$CD.m==max(W$CD.m))]
    STAT.m = max(W$CD.m)
    COMP   = W$Comparison[which(W$CD.m==max(W$CD.m))]
    names(STAT)   = "CD"
    names(STAT.m) = "CD.m"
    names(COMP)   = "Comparison"
  }
  
  if(statistic == "r"){
    STAT   = W$r[which(W$r.m==max(W$r.m))]
    STAT.m = max(W$r.m)
    COMP   = W$Comparison[which(W$r.m==max(W$r.m))]
    names(STAT)   = "r"
    names(STAT.m) = "r.m"
    names(COMP)   = "Comparison"
  }
  
  Z = list(pairs = W,
           comparison = COMP,
           statistic = STAT,
           statistic.m = STAT.m)
    
  return(Z)
  }
