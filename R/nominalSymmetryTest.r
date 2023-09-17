#' @title Exact and McNemar symmetry tests for paired contingency tables
#'
#' @description Conducts an omnibus symmetry test for a paired 
#'              contingency table and then post-hoc pairwise tests.  
#'              This is similar to  
#'              McNemar and McNemar-Bowker tests in use.
#' 
#' @param x A two-way contingency table. It must be square. 
#'          It can have two or
#'          more levels for each dimension.
#' @param method The method to adjust multiple p-values. 
#'               See \code{stats::p.adjust}.
#' @param exact If \code{TRUE}, uses the \code{binom.test} function.
#'              If \code{FALSE}, uses the \code{mcnemar.test} function.               
#' @param digits The number of significant digits in the output.
#' @param ... Additional arguments
#' 
#' @details The omnibus McNemar test
#'          may fail when there are zeros in critical
#'          cells.
#'          
#'          Currently, the \code{exact=TRUE} with a table greater
#'          than 2 x 2 will not produce an omnibus test result. 
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{https://rcompanion.org/handbook/H_05.html}
#' 
#' @seealso \code{\link{pairwiseMcnemar}}, 
#'          \code{\link{groupwiseCMH}},
#'          \code{\link{pairwiseNominalIndependence}}, 
#'          \code{\link{pairwiseNominalMatrix}}
#' 
#' @concept post-hoc
#' @concept McNemar's test
#' @concept McNemar Bowker
#' 
#' @return A list containing: a data frame of results of the global test;
#'         a data frame of results of the pairwise results;
#'         and a data frame mentioning the p-value adjustment method.
#'         
#' @examples
#' ### 2 x 2 repeated matrix example
#' data(AndersonRainBarrel)
#' nominalSymmetryTest(AndersonRainBarrel)
#'                     
#' ### 3 x 3 repeated matrix example
#' data(AndersonRainGarden)
#' nominalSymmetryTest(AndersonRainGarden,
#'                     exact = FALSE)
#'                                                               
#' @importFrom stats binom.test p.adjust
#' 
#' @export

nominalSymmetryTest = 
  function(x, method="fdr", digits=3, exact=FALSE, ...)
  {
     n = nrow(x)
     m = ncol(x)
     if((n < 2) | (m != n)){
        stop("Matrix must be square with at least two rows and columns")}
     N = n*n-n
     Y = rep(1/N, N)
     X = rep(NA, N)
     k=0
     for(i in 1:(n)){
        for(j in 1:n){
        if(i != j){
           k=k+1
           X[k]=x[i,j]
          }
         }
     }
  if(n>2){
    if(exact){
      p.value=NA
      }
    if(!exact){
      p.value=signif(mcnemar.test(x)$p.value, digits=digits)
      }
     cat("\n")
     }
  if(n==2){
     if(exact){
       p.value=signif(binom.test(X[1], sum(X), 0.5)$p.value, digits=digits)
       }
    if(!exact){
      p.value=signif(mcnemar.test(x)$p.value, digits=digits)
      }
     cat("\n")}
  Dimensions = paste(n, "x", n)
  if(n==2){
     V = data.frame(Dimensions, p.value)
     if(exact){Method="binomial test"}
     if(!exact){Method="McNemar test"}
     S = data.frame(Method)
     W = list(V, S)
     names(W) = c("Global.test.for.symmetry",
                  "Statistical.method")
     return(W)
   }   
  if(n>2){
     N = n*(n-1)/2
     Z = data.frame(Comparison=rep("A", N),
                    p.value=rep(NA, N),
                    p.adjust=rep(NA, N),
                    stringsAsFactors=FALSE)
     k=0               
     for(i in 1:(n-1)){
        for(j in (i+1):n){
           k=k+1
           Namea = as.character(rownames(x)[i])
           Nameb = as.character(colnames(x)[i])
           Namec = as.character(rownames(x)[j])
           Named = as.character(colnames(x)[j])
           a = x[i,j]
           b = x[j,i]
           d = x[i,j]+x[j,i]
           if(d>0){
             if(exact){
               P = signif(binom.test(a, d, 0.5)$p.value, digits=digits)
               }
             if(!exact){
               M=matrix(c(0, a, b, 0), nrow=2)
               P = signif(mcnemar.test(M)$p.value, digits=digits)
               }
             }
           if(d==0){P = NA}
           if(d<0){P = NA}
           P.adjust = NA
           Z[k,] =c( paste0(Namea, "/", Nameb, " : ", Namec, "/", Named), 
                     P, P.adjust)
        }
     }
     Z$p.adjust = signif(p.adjust(Z$p.value, method = method), digits=digits) 
     V = data.frame(Dimensions, p.value)
     U = data.frame(Method=method)
     if(exact){Method="binomial test"}
     if(!exact){Method="McNemar test"}
     S = data.frame(Method)
     W = list(V, Z, U, S)
     names(W) = c("Global.test.for.symmetry",
                  "Pairwise.symmetry.tests",
                  "p.adjustment",
                  "statistical.method")
     return(W)
  }
  }
