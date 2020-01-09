#' @title Matched-pairs rank biserial correlation coefficient
#' 
#' @description Calculates matched-pairs rank biserial correlation coefficient
#'              effect size
#'              for paired Wilcoxon signed-rank test; confidence intervals
#'              by bootstrap.
#' 
#' @param x A vector of observations.
#' @param g The vector of observations for
#'          the grouping, nominal variable.
#'          Only the first two levels of the nominal variable are used.
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
#' @param digits The number of significant digits in the output.
#' @param verbose If \code{TRUE}, prints information on sample size and ranks.
#' @param ... Additional arguments passed to \code{rank}             
#' 
#' @details  It is recommended that \code{NA}s be removed
#'           beforehand.
#'           
#'           When the data in the first group are greater than
#'           in the second group, rc is positive.
#'           When the data in the second group are greater than
#'           in the first group, rc is negative.
#'           Be cautious with this interpretation, as R will alphabetize
#'           groups if \code{g} is not already a factor.
#'           
#'           When rc is close to extremes,
#'           or with small counts in some cells,
#'           the confidence intervals 
#'           determined by this
#'           method may not be reliable, or the procedure may fail.
#'                      
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references King, B.M., P.J. Rosopa, and E.W. Minium. 2011. 
#'             Statistical Reasoning in the Behavioral Sciences, 6th ed.
#' @seealso \code{\link{wilcoxonPairedR}}
#' @concept correlation effect size ordinal nominal
#' @return A single statistic, rc.  
#'         Or a small data frame consisting of rc,
#'         and the lower and upper confidence limits.  
#'         
#' @examples
#' data(Pooh)
#' wilcox.test(Likert ~ Time, data=Pooh, paired=TRUE, exact=FALSE)
#' wilcoxonPairedRC(x = Pooh$Likert, g = Pooh$Time)
#' 
#' ### Example from King, Rosopa, and Minium
#' Placebo = c(24,39,29,28,25,32,31,33,31,22)
#' Drug    = c(28,29,34,21,28,15,17,28,16,12)
#' Y = c(Placebo, Drug)
#' Group = factor(c(rep("Placebo", length(Placebo)),  
#'                  rep("Drug", length(Drug))), 
#'                  levels=c("Placebo", "Drug"))
#' wilcoxonPairedRC(x = Y, g = Group)
#' 
#' ### Example with some zero differences
#' A = c(11,12,13,14,15,16,17,18,19,20)
#' B = c(12,14,16,18,20,22,12,10,19,20)
#' Y = c(A, B)
#' Group = factor(c(rep("A", length(A)),  
#'                  rep("B", length(B))))
#' wilcoxonPairedRC(x = Y, g = Group, verbose=TRUE, zero.method="Wilcoxon")
#' wilcoxonPairedRC(x = Y, g = Group, verbose=TRUE, zero.method="Pratt")
#' wilcoxonPairedRC(x = Y, g = Group, verbose=TRUE, zero.method="none")
#' 
#' @importFrom boot boot boot.ci
#' 
#' @export
 
wilcoxonPairedRC = function (x, g=NULL,
                      zero.method="Wilcoxon",
                      ci=FALSE, conf=0.95, type="perc",
                      R=1000, histogram=FALSE, digits=3,
                      verbose=FALSE,  ... ){
  
  if(is.factor(g)==F){g=factor(g)}
  x = x[as.numeric(g)<3]
  g = g[as.numeric(g)<3]
  g = droplevels(g)
  
  X1 = x[as.numeric(g)==1]
  X2 = x[as.numeric(g)==2]
  
  if(length(X1)!=length(X2)){stop("Groups not the same length", call.=FALSE)}
  Diff   = X1 - X2
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
  cat("Levels:",levels(g))
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
  Data = data.frame(x1=x[as.numeric(g)==1], x2=x[as.numeric(g)==2])
  Function = function(input, index){
    Input = input[index,]
    Diff   = Input$x1 - Input$x2
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
  Boot = boot(Data, Function, R=R)
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