#' @title Post-hoc tests for Cochran-Mantel-Haenszel test
#'
#' @description Conducts groupwise tests of association on a three-way 
#'              contingency table.
#' 
#' @param x A three-way contingency table.
#' @param group The dimension of the table to use as the grouping variable.
#'              Will be \code{1}, \code{2}, or \code{3}.
#' @param fisher If \code{TRUE}, conducts Fisher exact test.
#' @param gtest  If \code{TRUE}, conducts G test of association.
#' @param chisq  If \code{TRUE}, conducts Chi-square test of association.
#' @param method The method to use to adjust p-values.  See \code{?p.adjust}.
#' @param correct The correction to apply to the G test. 
#'                See \code{GTest}.
#' @param digits The number of digits for numbers in the output.
#' @param ... Other arguments passed to \code{chisq.test} or \code{GTest}.
#' 
#' @details If more than one of \code{fisher}, \code{gtest}, or \code{chisq} is
#'          set to \code{TRUE}, only one type of test of association
#'          will be conducted.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/H_06.html}
#' @concept Cochran Mantel Haenszel post-hoc fisher Chi square G
#' @return A data frame of groups, test used, p-values, and adjusted p-values.
#' 
#' @seealso \code{\link{nominalSymmetryTest}}, \code{\link{pairwiseMcnemar}},
#'          \code{\link{pairwiseNominalIndependence}}, 
#'          \code{\link{pairwiseNominalMatrix}} 
#'         
#' @examples
#' ### Post-hoc for Cochran-Mantel-Haenszel test
#' data(AndersonBias)
#' Table = xtabs(Count ~ Sex + Result + County,
#'               data=AndersonBias)
#' ftable(Table)
#' mantelhaen.test(Table)
#' groupwiseCMH(Table,
#'              group   = 3,
#'              fisher  = TRUE,
#'              gtest   = FALSE,
#'              chisq   = FALSE,
#'              method  = "fdr",
#'              correct = "none",
#'              digits  = 3)
#'                       
#' @importFrom stats mantelhaen.test chisq.test p.adjust
#' @importFrom DescTools GTest
#' 
#' @export

groupwiseCMH = 
  function(x, group=3,
           fisher=TRUE, gtest=FALSE, chisq=FALSE,
           method="fdr", correct="none", digits=3, ...) 
  {
     N = as.numeric(dim(x)[group])
     G = dimnames(x)[[group]]
     Y = data.frame(Group   = rep("A",N),
                    Test    = rep("A",N),
                    p.value = rep(NA,N),
                    adj.p   = rep(NA,N),
                    stringsAsFactors = FALSE)
     for(i in 1:N){
        if(group==1){X = x[i,,]}
        if(group==2){X = x[,i,]}
        if(group==3){X = x[,,i]}
     if(fisher==TRUE){
        Y[i,1] = G[i]
        Y[i,2] = "Fisher"
        Y[i,3] = signif(fisher.test(X, ...)$p.value, digits=digits)
        }
     if(gtest==TRUE){
        Y[i,1] = G[i]
        Y[i,2] = "G.test"
        Y[i,3] = signif(GTest(X, correct=correct, ...)$p.value, digits=digits)
        }
     if(chisq==TRUE){
        Y[i,1] = G[i]
        Y[i,2] = "Chi.square"
        Y[i,3] = signif(chisq.test(X, ...)$p.value, digits=digits)
        }
     }
     Y[,4] =  signif(p.adjust(Y[,3], method=method), digits=digits)
cat("\n")     
return(Y)   
}
