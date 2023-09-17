#' @title Pairwise tests of independence for nominal data with matrix output
#'
#' @description Conducts pairwise tests for a 2-dimensional matrix,
#'              in which at at least one dimension has more than two
#'              levels, as a post-hoc test. Conducts Fisher exact, Chi-square,
#'              or G-test.
#' 
#' @param x A two-way contingency table. At least one dimension should have
#'          more than two levels.
#' @param compare If \code{"row"}, treats the rows as the grouping variable.
#'                If \code{"column"}, treats the columns as the grouping variable.
#' @param fisher  If \code{"TRUE"}, conducts fisher exact test.
#' @param gtest   If \code{"TRUE"}, conducts G-test.
#' @param chisq   If \code{"TRUE"}, conducts Chi-square test of association.
#' @param method  The method to adjust multiple p-values. 
#'                See \code{\link{p.adjust}}.
#' @param correct The correction method to pass to \code{DescTools::GTest}.
#' @param digits The number of significant digits in the output.
#' @param ... Additional arguments, passed to \code{stats::fisher.test}, 
#'            \code{DescTools::GTest}, or \code{stats::chisq.test}.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{https://rcompanion.org/handbook/H_04.html}
#' 
#' @seealso \code{\link{pairwiseMcnemar}}, 
#'          \code{\link{groupwiseCMH}},
#'          \code{\link{nominalSymmetryTest}}, 
#'          \code{\link{pairwiseNominalIndependence}}
#'            
#' @concept post-hoc
#' @concept chi square test
#' @concept G test
#' @concept Fisher's exact test
#' 
#' @return A list consisting of: the test used,
#'         a matrix of unadjusted p-values,
#'         the p-value adjustment method used,
#'         and a matrix of adjusted p-values.
#'         
#' @examples
#' ### Independence test for a 4 x 2 matrix
#' data(Anderson)
#' fisher.test(Anderson)
#' Anderson = Anderson[(c("Heimlich", "Bloom", "Dougal", "Cobblestone")),]
#' PT = pairwiseNominalMatrix(Anderson,
#'                            fisher = TRUE,
#'                            gtest  = FALSE,
#'                            chisq  = FALSE)$Adjusted
#' PT
#' library(multcompView)
#' multcompLetters(PT)
#'                                                               
#' @importFrom stats fisher.test chisq.test
#' @importFrom DescTools GTest
#' @importFrom multcompView multcompLetters
#' 
#' @export

pairwiseNominalMatrix = 
  function(x, compare="row",
           fisher=TRUE, gtest=FALSE, chisq=FALSE,
           method="fdr", correct="none", digits=3, ...) 
  {
  if(compare=="row"){n = nrow(x)}
  if(compare=="column" | compare=="col"){n = ncol(x)}
     N = n*n
     Y = matrix(rep(NA_real_, N),ncol=n)
     Z = matrix(rep(NA_real_, N),ncol=n)
     if (compare=="row"){  
        rownames(Y)=rownames(x)
        colnames(Y)=rownames(x)
        rownames(Z)=rownames(x)
        colnames(Z)=rownames(x)
        }
     if (compare=="col" | compare=="column"){  
        rownames(Y)=colnames(x)
        colnames(Y)=colnames(x)
        rownames(Z)=colnames(x)
        colnames(Z)=colnames(x)
        }  
     k=0
     for(i in 1:(n-1)){
        for(j in (i+1):n){
           k=k+1
     Dataz = matrix(c(x[i,],x[j,]), nrow=2, byrow=TRUE)
     if(fisher==TRUE){
        Y[i,j] = signif(fisher.test(Dataz, ...)$p.value, digits=digits)}
     if(gtest==TRUE){
        Y[i,j] = signif(GTest(Dataz, correct=correct, ...)$p.value, digits=digits)}
     if(chisq==TRUE){
        Y[i,j] = signif(chisq.test(Dataz, ...)$p.value, digits=digits)}       
     }
    }
 if(fisher==TRUE){X = "Fisher exact test"}
 if(gtest ==TRUE){X = "G-test"}
 if(chisq ==TRUE){X = "Chi-square test"}
 Z[upper.tri(Z)] = 
        signif(p.adjust(Y[upper.tri(Y)], method=method), digits=digits)
 Z = t(Z)
 Z[upper.tri(Z)] = 
        signif(p.adjust(Y[upper.tri(Y)], method=method), digits=digits)
 diag(Z) = signif(1.00, digits = digits)

 W = method
 V = list(X, Y, W, Z)
 
     names(V) = c("Test",
                  "Unadjusted",
                  "Method",
                  "Adjusted")
return(V)   
} 
