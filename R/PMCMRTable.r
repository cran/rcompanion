#' @title Convert PMCMR Objects to a Data Frame
#'
#' @description Extracts a data frame of comparisons and p-values
#'              from an PMCMR object from the PMCMRplus
#'              package
#'
#' @param PMCMR   A PMCMR object
#' @param reverse If \code{TRUE}, reports the comparison as e.g. (B - A = 0).
#'                This will more closely match the output of 
#'                \code{PMCMRplus::summary.PMCMR} for all-pairs comparisons.
#'                If \code{FALSE}, reports the comparison as e.g. (A - B = 0).
#'                This will result in the output from \code{rcompanion::cldList}
#'                matching the output of
#'                \code{PMCMRplus::summaryGroup}
#'                              
#' @param digits  The significant digits in the output
#'
#' @details Should produce meaningful output for all-pairs 
#'          and many-to-one comparisons.          
#'
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_08.html}
#' @return A data frame of comparisons and p-values
#' 
#' @importFrom stats na.omit
#' 
#' @export

PMCMRTable = function(PMCMR, reverse=TRUE, digits=3){

G = PMCMR$model$g
P = signif(na.omit(as.vector(PMCMR$p.value)), digits=digits)

n = length(levels(G))
  N = n*(n-1)/2
  Z = data.frame(Comparison=rep("A", N),
                 p.value=rep(NA, N),
                 stringsAsFactors=FALSE)
  k=0               
  for(i in 1:(n-1)){
     for(j in (i+1):n){
       k=k+1
     Namea = as.character(levels(G)[i])
     Nameb = as.character(levels(G)[j])
     if(reverse==TRUE){ Z[k,] =c( paste0(Nameb, " - ", Namea, " = 0"), P[k])}
     if(reverse==FALSE){Z[k,] =c( paste0(Namea, " - ", Nameb, " = 0"), P[k])}
     }
  }
  return(na.omit(Z))
}