#' @title Rank biserial correlation coefficient for one-sample Wilcoxon test
#' 
#' @description Calculates rank biserial correlation coefficient
#'              effect size
#'              for one-sample Wilcoxon signed-rank test; confidence intervals
#'              by bootstrap.
#' 
#' @param x A vector of observations.
#' @param mu The value to compare \code{x} to, as in \code{wilcox.test}
#' @param zero.method If \code{"Wilcoxon"},
#'                    differences of zero are discarded and then ranks
#'                    are determined.
#'                    If \code{"Pratt"},
#'                    ranks are determined, 
#'                    and then differences of zero are discarded.
#'                     If \code{"none"},
#'                    differences of zero are not discarded.
#' @param ci If \code{TRUE}, returns confidence intervals by bootstrap.
#'           May be slow.
#' @param conf The level for the confidence interval.
#' @param type The type of confidence interval to use.
#'             Can be any of "\code{norm}", "\code{basic}", 
#'                           "\code{perc}", or "\code{bca}".
#'             Passed to \code{boot.ci}.
#' @param R The number of replications to use for bootstrap.
#' @param histogram If \code{TRUE}, produces a histogram of bootstrapped values.
#' @param verbose If \code{TRUE}, prints information on sample size and ranks.
#' @param digits The number of significant digits in the output.
#' @param ... Additional arguments passed to the \code{wilcoxsign_test} function.
#' 
#' @details  It is recommended that \code{NA}s be removed
#'           beforehand.
#'           
#'           When rc is close to extremes,
#'           or with small counts in some cells,
#'           the confidence intervals 
#'           determined by this
#'           method may not be reliable, or the procedure may fail.
#'                      
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_02.html}
#' @seealso \code{\link{wilcoxonPairedRC}}
#' @concept correlation effect size ordinal nominal
#' @return A single statistic, rc.  
#'         Or a small data frame consisting of rc,
#'         and the lower and upper confidence limits.  
#'         
#' @examples
#' ### Example with one zero difference
#' A = c(11,12,13,14,15,16,17,18,19,20)
#' #' wilcoxonOneSampleRC(x = A, mu=15)
#' wilcoxonOneSampleRC(x = A, mu=15, verbose=TRUE, zero.method="Wilcoxon")
#' wilcoxonOneSampleRC(x = A, mu=15, verbose=TRUE, zero.method="Pratt")
#' wilcoxonOneSampleRC(x = A, mu=15, verbose=TRUE, zero.method="none")
#' 
#' @importFrom boot boot boot.ci
#' 
#' @export
 
wilcoxonOneSampleRC = function (x, mu=NULL,
                      zero.method="Wilcoxon",
                      ci=FALSE, conf=0.95, type="perc",
                      R=1000, histogram=FALSE, digits=3, 
                      verbose=FALSE,  ... ){
  
  Diff   = x - mu
  Diff2  = abs(Diff)

  if(zero.method=="none")
    {RR     = -1 * rank(Diff2) * sign(Diff)}
  if(zero.method=="Pratt")
    {R      = -1 * rank(Diff2) * sign(Diff)
     RR     = R[R!=0]}
  if(zero.method=="Wilcoxon")
    {Diff   = Diff[Diff!=0]
     Diff2  = Diff2[Diff2!=0]
     RR     = -1 * rank(Diff2) * sign(Diff)
     }

  Rplus  = sum(RR[RR>0])
  Rminus = sum(abs(RR[RR<0]))
  Tee    = min(Rplus, Rminus)
  n      = length(RR)
  if(Rplus>=Rminus){
             RC = -4 * abs((Tee - (Rplus + Rminus)/2) / n / (n+1))}
  if(Rplus<Rminus){
             RC = 4 * abs((Tee - (Rplus + Rminus)/2) / n / (n+1))}
  RC     = signif(RC, digits=digits)
  
  if(verbose)
  {
  cat("\n")
  cat("zero.method:", zero.method)
  cat("\n")
  cat("n kept", "=", length(RR))
  cat("\n")
  cat("Ranks plus", "=", Rplus)
  cat("\n")
  cat("Ranks minus", "=", Rminus)
  cat("\n")
  cat("T value =",  Tee)
  cat("\n")
  cat("\n")
  }
    
if(ci==TRUE){
  Function = function(input, index){
    Input = input[index]
    Diff   = Input - mu
    Diff2  = abs(Diff)
    if(zero.method=="none")
      {RR     = -1 * rank(Diff2) * sign(Diff)}
    if(zero.method=="Pratt")
      {R      = -1 * rank(Diff2) * sign(Diff)
       RR     = R[R!=0]}
    if(zero.method=="Wilcoxon")
      {Diff   = Diff[Diff!=0]
       Diff2  = Diff2[Diff2!=0]
       RR     = -1 * rank(Diff2) * sign(Diff)
       }
    Rplus  = sum(RR[RR>0])
    Rminus = sum(abs(RR[RR<0]))
    Tee    = min(Rplus, Rminus)
    n      = length(RR)
    if(Rplus>=Rminus){
               RC = -4 * abs((Tee - (Rplus + Rminus)/2) / n / (n+1))}
    if(Rplus<Rminus){
               RC = 4 * abs((Tee - (Rplus + Rminus)/2) / n / (n+1))}
    return(RC)                  
  }
  Boot = boot(x, Function, R=R)
  BCI  = boot.ci(Boot, conf=conf, type=type)
  if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3]}
  if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5]}
  if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5]}
  if(type=="bca")  {CI1=BCI$bca[4];     CI2=BCI$bca[5]}  
  
  CI1=signif(CI1, digits=digits)
  CI2=signif(CI2, digits=digits)
  
  if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray", xlab="rc",
                                       main="")}
  }
if(ci==FALSE){names(RC)="rc"; return(RC)}
if(ci==TRUE){DF=data.frame(rc=RC, lower.ci=CI1, upper.ci=CI2)
             rownames(DF) = 1:nrow(DF)
             return(DF)
             }  
}