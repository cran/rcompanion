#' @title Cohen's h to compare proportions for 2 x 2 contingency tables
#'
#' @description Calculates Cohen's h
#'              for 2 x 2 contingency tables, such as those that 
#'              might be analyzed with   
#'              a chi-square test of association.
#' 
#' @param x A 2 x 2 contingency table.
#' @param observation If \code{"row"}, the row constitutes an observation.  
#'                    That is, the sum of each row is 100 percent.
#'                    If \code{"column"}, the column constitutes an observation.  
#'                    That is, the sum of each column is 100 percent.
#' @param verbose If \code{TRUE}, prints the proportions for each observation.
#' @param digits The number of significant digits in the output.
#' 
#' @details Cohen's h is an effect size to compare two proportions.
#'          For a 2 x 2 table: 
#'          Cohen's h equals Phi2 - Phi1, where,
#'          If observations are in rows, P1 = a/(a+b) and P2 = c/(c+d).
#'          If observations are in columns, P1 = a/(a+c) and P2 = b/(b+d).          
#'          Phi = 2 * asin(sqrt(P))
#.
#'        
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/H_05.html}
#' @concept Effect size Cohen h contingency table nominal
#' @return A single statistic.
#'         
#' @seealso \code{\link{cohenG}}
#'         
#' @examples
#' data(Pennsylvania18)
#' Pennsylvania18
#' cohenH(Pennsylvania18, observation="row")
#' 
#' @export

cohenH = 
  function(x, observation="row", verbose = TRUE, digits=3)
  {
     cat("\n")
     n = nrow(x)
     m = ncol(x)
     if((n != 2) | (m != n)){
        stop("Matrix must be square with exactly two rows and columns")}
     a = x[1,1]
     b = x[1,2]
     c = x[2,1]
     d = x[2,2]
     if(observation=="row"){
       P1 = a/(a+b)
       P2 = c/(c+d)
       if(verbose){
         if(is.null(rownames(x)[1]) | is.null(rownames(x)[2])){
           rownames(x)=c("Row 1", "Row 2")}
         Z = data.frame(Group = c(rownames(x)[1], rownames(x)[2]),
                        Proportion = 
                          c(signif(P1,digits=digits), signif(P2,digits=digits)))
         print(Z)
         cat("\n")
       }
     }
     if(observation=="column"){
       P1 = a/(a+c)
       P2 = b/(b+d)
       if(verbose){
         if(is.null(colnames(x)[1]) | is.null(colnames(x)[2])){
           colnames(x)=c("Column 1", "Column 2")}
         Z = data.frame(Group = c(colnames(x)[1], colnames(x)[2]),
                        Proportion = 
                          c(signif(P1,digits=digits), signif(P2,digits=digits)))
         print(Z)
         cat("\n")
       }
     }
     Phi1 = 2 * asin(sqrt(P1))
     Phi2 = 2 * asin(sqrt(P2))
     Cohen.h = signif((Phi2 - Phi1), digits=digits)
     names(Cohen.h) = c("Cohen's h")
     return(Cohen.h)
  }
  
 