#' @title Pairwise two-sample independence tests with matrix output
#'
#' @description Conducts pairwise two-sample independence tests across groups.
#' 
#' @param formula A formula indicating the measurement variable and
#'                the grouping variable. e.g. y ~ group.
#' @param data   The data frame to use.                
#' @param x      The response variable as a vector.
#' @param g      The grouping variable as a vector.
#' @param method The p-value adjustment method to use for multiple tests.
#'               See \code{stats::p.adjust}.
#' @param ...    Additional arguments passed to
#'               \code{coin::independence_test}.               
#'             
#' @details The input should include either \code{formula} and \code{data};
#'          or \code{x}, and \code{g}.
#' 
#'          This function is a wrapper for \code{coin::independence_test},
#'          passing pairwise groups to the function. It's critical to read
#'          and understand the documentation for this function to understand
#'          its use and options.
#'          
#'          For some options for common tests, see Horthorn et al., 2008.
#'          
#' @note    The parsing of the formula is simplistic. 
#'          The first variable on the
#'          left side is used as the measurement variable.  
#'          The first variable on the
#'          right side is used for the grouping variable.                                                                                 
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{http://rcompanion.org/handbook/K_02.html}
#' 
#'             Hothorn, T., K. Hornik, M.A. van de Wiel, and A. Zeileis. 2008. 
#'             Implementing a Class of Permutation Tests: The coin Package. 
#'             Journal of Statistical Software, 28(8), 1–23.
#'             
#' @seealso \code{\link{pairwisePermutationTest}}
#' 
#' @concept post-hoc 
#' @concept permutation test
#' 
#' @return A list consisting of:
#'         A matrix of p-values;
#'         the p-value adjustment method;
#'         a matrix of adjusted p-values.
#'
#' @examples
#' ### Fisher-Pitman test
#' 
#' data(BrendonSmall)
#' 
#' library(coin)
#'                                  
#' independence_test(Sodium ~ Instructor, data = BrendonSmall, 
#'                   teststat = "quadratic") 
#'                                       
#' PT = pairwisePermutationMatrix(Sodium ~ Instructor,
#'                                data     = BrendonSmall,
#'                                teststat = "quadratic",
#'                                method   = "fdr")
#' PT
#' 
#' PA = PT$Adjusted
#' library(multcompView)
#' multcompLetters(PA,
#'                 compare="<",
#'                 threshold=0.05,
#'                 Letters=letters)   
#' 
#' @importFrom stats p.adjust
#' @importFrom coin independence_test statistic
#' 
#' @export  

pairwisePermutationMatrix = 
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
     z = independence_test(x ~ g2, data=Dataz, ...)                  
   Y[i,j] = signif(pvalue(z), digits = 4)
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