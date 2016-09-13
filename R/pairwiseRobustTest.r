#' @title Pairwise two-sample robust tests
#'
#' @description Performs pairwise two-sample robust tests across groups.
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
#'          The \code{pairwiseRobustTest} function 
#'          can be used as a post-hoc method following an omnibus 
#'          one-way anova with robust estimation.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/rcompanion/d_08a.html}
#' @seealso \code{\link{pairwiseRobustMatrix}}
#' @concept robust Huber post-hoc one-way
#' @return A dataframe of the groups being compared, the p-values,
#'         and the adjusted p-values. 
#'         
#' @examples
#' data(PoohPiglet)
#' PoohPiglet = PoohPiglet[order(factor(PoohPiglet$Speaker, 
#'                         levels=c("Pooh", "Tigger", "Piglet"))),]
#' pairwiseRobustTest(x      = PoohPiglet$Likert,
#'                    g      = PoohPiglet$Speaker,
#'                    method = "fdr")
#' 
#' @importFrom stats p.adjust anova
#' @importFrom WRS2 pb2gen
#' 
#' @export

pairwiseRobustTest = 
  function(x, g, est="mom", nboot=599, method="fdr", ...)
  {
  n = length(unique(g))
  N = n*(n-1)/2
  d = data.frame(x = x, g = g)
  Z = data.frame(Comparison=rep("A", N),
                 Statistic=rep(NA, N),
                 p.value=rep(NA, N),
                 p.adjust=rep(NA, N),
                 stringsAsFactors=FALSE)
  k=0
  for(i in 1:(n-1)){
     for(j in (i+1):n){
       k=k+1
     Namea = as.character(unique(g)[i])
     Nameb = as.character(unique(g)[j])
     Datax = subset(d, g==unique(g)[i])
     Datay = subset(d, g==unique(g)[j])
     Dataz = rbind(Datax, Datay)
     Dataz$g2 = factor(Dataz$g)
     print(paste0("comparison ",k," ..."))
     z = pb2gen(x ~ g2, data=Dataz,
                      est=est,
                      nboot=nboot, 
                      ...)                           
     P = signif(z$p.value, digits=4)
     S = signif(z$test, digits=4)
     P.adjust = NA    
     Z[k,] =c( paste0(Namea, " - ", Nameb, " = 0"), 
             S, P, P.adjust)           
     }
    } 
  Z$p.adjust = 
      signif(p.adjust(Z$p.value, method = method), digits=4) 
  cat("\n","\n")
  return(Z)
  }