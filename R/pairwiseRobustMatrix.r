#' @title Pairwise two-sample robust tests with matrix output
#'
#' @description Performs pairwise two-sample robust tests across groups
#'              with matrix output.
#' 
#' @param x      The response variable as a vector.
#' @param g      The grouping variable as a vector.
#' @param est    Estimate used for group comparisons.
#'               \code{"onestep"}, \code{"mom"}, \code{"median"}, 
#'               or \code{"mean"}.
#'               See \code{\link{pb2gen}} for details.
#' @param nboot  The number of bootstrap samples.
#' @param method The p-value adjustment method to use for multiple tests.
#'               See \code{\link{p.adjust}}.
#' @param ...    Additional arguments passed to
#'               \code{\link{pb2gen}}.               
#'             
#' @details 
#'          The \code{WRS2} package provides functions for robust estimation
#'          and hypothesis testing.  This function invokes the
#'          \code{pb2gen} to make pairwise comparisons among
#'          groups.
#'          
#'          The \code{pairwiseRobustMatrix} function 
#'          can be used as a post-hoc method following an omnibus 
#'          one-way anova with robust estimation.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/rcompanion/d_08a.html}
#' @seealso \code{\link{pairwiseRobustTest}}
#' @concept robust Huber post-hoc one-way
#' @return A list consisting of:
#'         A matrix of p-values;
#'         the p-value adjustment method;
#'         a matrix of adjusted p-values. 
#'         
#' @examples
#' data(PoohPiglet)
#' PoohPiglet$Speaker = factor(PoohPiglet$Speaker, 
#'                             levels=c("Pooh", "Tigger", "Piglet"))               
#' PT = pairwiseRobustMatrix(x      = PoohPiglet$Likert,
#'                           g      = PoohPiglet$Speaker,
#'                           method = "fdr")$Adjusted
#' PT                           
#' library(multcompView)
#' multcompLetters(PT,
#'                 compare="<",
#'                 threshold=0.05,
#'                 Letters=letters)
#' 
#' @importFrom stats p.adjust anova
#' @importFrom WRS2 pb2gen
#' 
#' @export

pairwiseRobustMatrix = 
  function(x, g, est="mom", nboot=599, method="fdr", ...)
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
     cat("comparison ",k," ...","\n")
     Datax = subset(d, g==levels(g)[i])
     Datay = subset(d, g==levels(g)[j])
     Dataz = rbind(Datax, Datay)
     Dataz$g2 = factor(Dataz$g)
     z = pb2gen(x ~ g2, data=Dataz,
                      est=est,
                      nboot=nboot,
                      ...)                   
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
