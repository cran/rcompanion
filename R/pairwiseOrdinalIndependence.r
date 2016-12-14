#' @title Pairwise tests of independence for tables with one ordinal variable
#'
#' @description Conducts pairwise tests for a 2-dimensional table,
#'              in which one variable is ordinal and one variable
#'              is nominal.  The function relies on the \code{coin} package.
#' 
#' @param x A two-way contingency table. One dimension is ordinal and one
#'                is un-ordered nominal.
#' @param compare If \code{"row"}, treats the rows as the grouping variable.
#'                If \code{"column"}, treats the columns as the grouping 
#'                variable.
#' @param scores  Optional vector to specify the spacing of the ordinal
#'                variable.
#' @param method  The method to adjust multiple p-values. 
#'                See \code{\link{p.adjust}}.
#' @param digits The number of significant digits in the output.
#' @param ... Additional arguments, passed to \code{\link{chisq_test}}. 
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/H_09.html}
#' @concept Linear-by-linear Contingency table ordinal Chi-square.
#' @return A data frame of comparisons, p-values, and adjusted p-values.
#'         
#' @seealso \code{\link{pairwiseNominalIndependence}}
#'         
#' @examples
#' ### Independence test for table with one ordinal variable
#' data(Breakfast)
#' require(coin)
#' chisq_test(Breakfast,
#'            scores = list("Breakfast" = c(-2, -1, 0, 1, 2)))
#' PT = pairwiseOrdinalIndependence(Breakfast, compare = "row")
#' PT
#' cldList(comparison = PT$Comparison, 
#'         p.value    = PT$p.value, 
#'         threshold  = 0.05)
#'         
#' ### Similar to Kruskal-Wallis test for Likert data
#' data(PoohPiglet)
#' XT = xtabs(~ Speaker + Likert, data = PoohPiglet)
#' XT
#' require(coin)
#' chisq_test(XT,
#'            scores = list("Likert" = c(1, 2, 3, 4, 5)))
#' PT=pairwiseOrdinalIndependence(XT, compare = "row")
#' PT
#' cldList(comparison = PT$Comparison, 
#'         p.value    = PT$p.value, 
#'         threshold  = 0.05)         
#'                                                               
#' @importFrom coin chisq_test
#' 
#' @export

pairwiseOrdinalIndependence = 
  function(x, compare="row", scores=NULL,
           method="fdr", digits=3, ...) 
  {
  if(compare=="row"){
      n    = as.numeric(nrow(x))
      Name = names(dimnames(x))[2]
      m    = as.numeric(ncol(x))
      }
  if(compare=="column" | compare=="col"){
      n    = as.numeric(ncol(x))
      Name = names(dimnames(x))[1]
      m    = as.numeric(nrow(x))
      }
  N = n*(n-1)/2
  Z = data.frame(Comparison=rep("A", N),
                 p.value=rep(NA, N),
                 p.adjust=rep(NA, N),
                 stringsAsFactors=FALSE)
  if(is.null(scores)){scores = 1:m-(m+1)/2}
  k=0               
  for(i in 1:(n-1)){
     for(j in (i+1):n){
       k=k+1
 if (compare=="row"){
     T1     = x[i,]
     T2     = x[j,]
     Matrix = rbind(T1, T2)
     Namea  = rownames(x)[i]
     Nameb  = rownames(x)[j]
     Table  = as.table(Matrix)
     names(dimnames(Table)) = c("Row", "Column")
     PV     = signif(pvalue(chisq_test(Table,
                                       scores=list(Column = scores),
                                       ...)),
                                       digits=digits)
     }
 if (compare=="column" | compare=="col"){  
     T1    = x[,i]
     T2    = x[,j]
     Matrix = cbind(T1, T2)
     Namea = colnames(x)[i]
     Nameb = colnames(x)[j]
     Table  = as.table(Matrix)
     names(dimnames(Table)) = c("Row", "Column")
     PV     = signif(pvalue(chisq_test(Table,
                                       scores=list(Row = scores),
                                       ...)),
                                       digits=digits)
     }     
 Z$Comparison[k] = paste0(Namea, " : ", Nameb)
 Z$p.value[k]    = PV  
 }
 }
 Z$p.adjust = signif(p.adjust(Z$p.value, method = method), digits=digits)
 
 return(Z)
}