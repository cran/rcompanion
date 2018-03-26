#' @title Cohen's g and odds ratio for paired contingency tables
#'
#' @description Calculates Cohen's g and odds ratio
#'              for paired contingency tables, such as those that 
#'              might be analyzed with   
#'              McNemar or McNemar-Bowker tests.
#' 
#' @param x A two-way contingency table. It must be square. 
#'          It can have two or
#'          more levels for each dimension.
#' @param digits The number of significant digits in the output.
#' 
#' @details For a 2 x 2 table, where a and d are the concordant cells
#'          and b and c are discordant cells:
#'          Odds ratio is the greater of b/c or c/b;
#'          P is the greater of b/(b+c) or c/(b+c);
#'          and Cohen's g is P - 0.5.
#'          These statistics are extended to tables larger
#'          than 2 x 2.
#'        
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/H_05.html}
#' @concept Effect size Cohen g McNemar Bowker contingency symmetry nominal
#' @return A list containing: a data frame of results of the global statistics;
#'         and a data frame of results of the pairwise statistics.
#'         
#' @seealso \code{\link{nominalSymmetryTest}}
#'         
#' @examples
#' ### 2 x 2 repeated matrix example
#' data(AndersonRainBarrel)
#' cohenG(AndersonRainBarrel)
#'                     
#' ### 3 x 3 repeated matrix
#' data(AndersonRainGarden)
#' cohenG(AndersonRainGarden)
#' 
#' @export

cohenG = 
  function(x, digits=3)
  {
     cat("\n")
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
    if(n==2){
      b  = x[1,2]
      c  = x[2,1]
      OR = max(b/c, c/b)
      OR = signif(OR, digits=digits)
      P  = max(b/(b+c), c/(c+b))
      P  = signif(P, digits=digits)
      g  = P - 0.5
      g  = signif(g, digits=digits)
      }

  Dimensions = paste(n, "x", n)
  
  if(n==2){
     V = data.frame(Dimensions, OR, P, g)
     W = list(V)
     names(W) = c("Global.statistics")
     return(W)
  }
  
  if(n>2){
     N = n*(n-1)/2
     Z = data.frame(Comparison = rep("A", N),
                    OR         = rep(NA, N),
                    P          = rep(NA, N),
                    g          = rep(NA, N),
                    stringsAsFactors=FALSE)
     k=0
     NUMERATOR = 0
     DENOMINATOR = 0
     DENOMINATOR1 = 0
     for(i in 1:(n-1)){
        for(j in (i+1):n){
           k=k+1
           Namea = as.character(rownames(x)[i])
           Nameb = as.character(colnames(x)[i])
           Namec = as.character(rownames(x)[j])
           Named = as.character(colnames(x)[j])
           b = x[i,j]
           c = x[j,i]
           OR = max(b/c, c/b)
           P  = max(b/(b+c), c/(c+b))
           g  = P - 0.5
           OR = signif(OR, digits=digits)
           P  = signif(P, digits=digits)
           g  = signif(g, digits=digits)
           Z[k,] =c( paste0(Namea, "/", Nameb, " : ", Namec, "/", Named),
                    OR, P, g)
           NUMERATOR   = NUMERATOR + max(b, c)
           DENOMINATOR = DENOMINATOR + min(b, c)
           DENOMINATOR1 = DENOMINATOR1 + b + c
        }
     }
     OR = NUMERATOR / DENOMINATOR
     P  = NUMERATOR / DENOMINATOR1
     g  = P - 0.5
     OR = signif(OR, digits=digits)
     P  = signif(P, digits=digits)
     g  = signif(g, digits=digits)
     
     V = data.frame(Dimensions, OR, P, g)
     W = list(V, Z)
     names(W) = c("Global.statistics",
                  "Pairwise.statistics")
     return(W)
  }
  }
