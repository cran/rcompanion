#' @title Pairwise tests of independence for nominal data
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
#' @param correct The correction method to pass to \code{GTest}
#'                See \code{\link{GTest}}.
#' @param digits The number of significant digits in the output.
#' @param ... Additional arguments, passed to \code{\link{fisher.test}}, 
#'            \code{\link{GTest}}, or \code{\link{chisq.test}}.
#' 
#' @details (None.)
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/H_04.html}
#' @concept Chi-square G-test Fisher contingency table nominal
#' @return A data frame of comparisons, p-values, and adjusted p-values.
#'         
#'  @seealso \code{\link{pairwiseMcnemar}}, \code{\link{groupwiseCMH}},
#'           \code{\link{nominalSymmetryTest}}, 
#'           \code{\link{pairwiseNominalMatrix}}
#'         
#' @examples
#' ### Independence test for a 4 x 2 matrix
#' data(Anderson)
#' fisher.test(Anderson)
#' pairwiseNominalIndependence(Anderson)
#'                                                               
#' @importFrom stats fisher.test chisq.test
#' @importFrom DescTools GTest
#' 
#' @export

pairwiseNominalIndependence = 
  function(x, compare="row",
           fisher=TRUE, gtest=TRUE, chisq=TRUE,
           method="fdr", correct="none", digits=3, ...) 
  {
  if(compare=="row"){n = nrow(x)}
  if(compare=="column" | compare=="col"){n = ncol(x)}
  N = n*(n-1)/2
  Z = data.frame(Comparison=rep("A", N),
                 stringsAsFactors=FALSE)
  p.Fisher = rep(NA, N)
  p.Gtest = rep(NA, N)
  p.Chisq = rep(NA, N)

  k=0               
  for(i in 1:(n-1)){
     for(j in (i+1):n){
       k=k+1
 if (compare=="row"){  
     Namea = as.character(rownames(x)[i])
     Nameb = as.character(rownames(x)[j])
     Dataz = matrix(c(x[i,],x[j,]), nrow=2, byrow=TRUE)
     }
 if (compare=="column" | compare=="col"){  
     Namea = as.character(colnames(x)[i])
     Nameb = as.character(colnames(x)[j])
     Dataz = matrix(c(x[,i],x[,j]), ncol=2, byrow=FALSE)
     }     
 Z$Comparison[k] = paste0(Namea, " : ", Nameb)     
 if(fisher==TRUE){
    p.Fisher[k] = signif(fisher.test(Dataz, ...)$p.value, digits=digits)}
 if(gtest==TRUE){
    p.Gtest[k] = signif(GTest(Dataz, correct=correct, ...)$p.value, digits=digits)}
 if(chisq==TRUE){
    p.Chisq[k] = signif(chisq.test(Dataz, ...)$p.value, digits=digits)}       
  }
 }
 if(fisher==TRUE){ 
    Z$p.Fisher = p.Fisher
    Z$p.adj.Fisher = signif(p.adjust(Z$p.Fisher, method = method), digits=digits)
  }
 if(gtest==TRUE){ 
    Z$p.Gtest = p.Gtest
    Z$p.adj.Gtest = signif(p.adjust(Z$p.Gtest, method = method), digits=digits)
  }
 if(chisq==TRUE){ 
     Z$p.Chisq = p.Chisq
     Z$p.adj.Chisq = signif(p.adjust(Z$p.Chisq, method = method), digits=digits)
  }
return(Z)
}