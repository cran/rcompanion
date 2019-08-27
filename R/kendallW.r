#' @title Kendall's W with bootstrapped confidence interval
#' 
#' @description Calculates Kendall's W coefficient of concordance,
#'              which can be used as an effect size statistic for
#'              unreplicated complete block design
#'              such as where Friedman's test might be used.
#'              This function is a wrapper for the \code{KendallW}
#'              function in the \code{DescTools} package,
#'              with the addition of bootstrapped
#'              confidence intervals.
#' 
#' @param x A k x m matrix or table, with k treatments in rows 
#'          and m raters or blocks in columns.
#' 
#' @param correct Passed to \code{KendallW}.
#' @param na.rm   Passed to \code{KendallW}.
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
#' @param ... Additional arguments passed to the \code{KendallW} function.             
#' 
#' @details  See the \code{KendallW} function in the \code{DescTools} package
#'           for details.
#'           
#'           When W is close to 0 or very large,
#'           or with small sample size,
#'           the confidence intervals 
#'           determined by this
#'           method may not be reliable, or the procedure may fail.
#'           
#'           Because W is always positive, if \code{type="perc"},
#'           the confidence interval will
#'           never cross zero, and should not
#'           be used for statistical inference.
#'           However, if \code{type="norm"}, the confidence interval
#'           may cross zero. 
#'           
#'           When producing confidence intervals by bootstrap, 
#'           this function treats each rater or block as an observation.
#'           It is not clear to the author if this approach produces accurate
#'           confidence intervals, but it appears to be reasonable.
#'           
#'                      
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_04.html}
#' @concept Friedman effect size
#' @return A single statistic, W.  
#'         Or a small data frame consisting of W,
#'         and the lower and upper confidence limits.  
#'         
#' @section Acknowledgments:
#'  Thanks to Indrajeet Patil, author of \code{ggstatsplot},
#'  and \code{groupedstats} for help in the inspiring and 
#'  coding of this function.
#'               
#' @examples
#' data(BobBelcher)
#' Table = xtabs(Likert ~ Instructor + Rater, data = BobBelcher)
#' kendallW(Table)
#' 
#' @importFrom DescTools KendallW
#' @importFrom boot boot boot.ci
#' 
#' @export
 
kendallW = function (x, correct=TRUE, na.rm=FALSE, 
                      ci=FALSE, conf=0.95, type="perc",
                      R=1000, histogram=FALSE, digits=3, ... ){

  KW = suppressWarnings(KendallW(x, correct=correct, na.rm=na.rm, ...))
  W  =signif(KW, digits=digits)

if(ci==TRUE){
  
  Function = function(input, index){
                    Input = input[index,]
                    KW = suppressWarnings(KendallW(t(Input), correct=correct, 
                                                   na.rm=na.rm, ...))
                    return(KW)}
  
  Boot = boot(t(x), Function, R=R)
  BCI  = boot.ci(Boot, conf=conf, type=type)
  if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3];}
  if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5];}
  if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5];}
  if(type=="bca")  {CI1=BCI$bca[4];     CI2=BCI$bca[5];}  
  
  CI1=signif(CI1, digits=digits)
  CI2=signif(CI2, digits=digits)
  
  if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray",
                      main="", xlab="W")}
}
  
if(ci==FALSE){names(W)="W"; return(W)}
if(ci==TRUE){return(data.frame(W=W, lower.ci=CI1, upper.ci=CI2))}  
  
}  
  
  