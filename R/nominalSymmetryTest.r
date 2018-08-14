#' @title Exact and Monte Carlo symmetry tests for paired contingency tables
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
#' @param digits The number of significant digits in the output.
#' @param ... Additional arguments, passed to \code{EMT::multinomial.test}.
#' 
#' @details If Monte Carlo is not used, the test of symmetry uses
#'          an exact test by conducting either a binomial 
#'          or multinomial goodness-of-fit test.
#'          
#'          These are equivalent to uncorrected 
#'          McNemar and McNemar-Bowker tests,
#'          but will not fail when there are zeros in critical
#'          cells, as will the \code{mcnemar.test} function.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/H_05.html}
#' @concept McNemar Bowker contingency table symmetry nominal
#' @return A list containing: a data frame of results of the global test;
#'         a data frame of results of the pairwise results;
#'         and a data frame mentioning the p-value adjustment method.
#'         
#' @seealso \code{\link{pairwiseMcnemar}}, \code{\link{groupwiseCMH}},
#'          \code{\link{pairwiseNominalIndependence}}, 
#'          \code{\link{pairwiseNominalMatrix}}
#'         
#' @examples
#' ### 2 x 2 repeated matrix example
#' data(AndersonRainBarrel)
#' nominalSymmetryTest(AndersonRainBarrel)
#'                     
#' ### 3 x 3 repeated matrix example with Monte Carlo
#' data(AndersonRainGarden)
#' nominalSymmetryTest(AndersonRainGarden,
#'                     MonteCarlo = TRUE,
#'                     ntrial     = 10000)
#'                     
#'### 4 x 4 repeated matrix example that fails with mcnemar.test
#' data(Religion)
#' nominalSymmetryTest(Religion,
#'                     MonteCarlo = TRUE,
#'                     ntrial     = 10000)
#'                                                               
#' @importFrom stats binom.test p.adjust
#' @importFrom EMT multinomial.test
#' 
#' @export

nominalSymmetryTest = 
  function(x, method="fdr", digits=3, ...)
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
     p.value=signif(multinomial.test(X,Y,...)$p.value, digits=digits)
     cat("\n")
     }
  if(n==2){
     p.value=signif(binom.test(X[1], sum(X), 0.5)$p.value, digits=digits)
     cat("\n")}
  Dimensions = paste(n, "x", n)
  if(n==2){
     V = data.frame(Dimensions, p.value)
     W = list(V)
     names(W) = c("Global.test.for.symmetry")
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
           b = x[i,j]+x[j,i]
           if(b>0){P = signif(binom.test(a, b, 0.5)$p.value, digits=digits)}
           if(b==0){P = NA}
           if(b<0){P = NA}
           P.adjust = NA
           Z[k,] =c( paste0(Namea, "/", Nameb, " : ", Namec, "/", Named), 
                     P, P.adjust)
        }
     }
     Z$p.adjust = signif(p.adjust(Z$p.value, method = method), digits=digits) 
     V = data.frame(Dimensions, p.value)
     U = data.frame(Method=method)
     W = list(V, Z, U)
     names(W) = c("Global.test.for.symmetry",
                  "Pairwise.symmetry.tests",
                  "p.adjustment")
     return(W)
  }
  }
