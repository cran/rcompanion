#' @title r effect size for Wilcoxon two-sample paired signed-rank test
#' 
#' @description Calculates r effect size
#'              for a Wilcoxon two-sample paired signed-rank test; 
#'              confidence intervals by bootstrap.
#' 
#' @param x A vector of observations.
#' @param g The vector of observations for
#'          the grouping, nominal variable.
#'          Only the first two levels of the nominal variable are used.
#'          The data must be ordered so that the first observation of the
#'          of the first group is paired with the first observation of the
#'          second group.
#' @param coin If \code{FALSE}, the default, the Z value
#'                is extracted from a function similar to the
#'                \code{wilcox.test} function in the stats package.
#'                If \code{TRUE}, the Z value
#'                is extracted from the \code{wilcox_test} function in the
#'                coin package.  This method may be much slower, especially
#'                if a confidence interval is produced.        
#' @param ci If \code{TRUE}, returns confidence intervals by bootstrap.
#'           May be slow.
#' @param conf The level for the confidence interval.
#' @param type The type of confidence interval to use.
#'             Can be any of "\code{norm}", "\code{basic}", 
#'                           "\code{perc}", or "\code{bca}".
#'             Passed to \code{boot.ci}.
#' @param R The number of replications to use for bootstrap.
#' @param histogram If \code{TRUE}, produces a histogram of bootstrapped values.
#' @param cases By default the \code{N} used in the formula for \code{r}
#'              is the number of pairs.  If \code{cases=FALSE},
#'              the \code{N} used in the formula for \code{r}
#'              is the total number of observations, as some sources suggest.
#' @param digits The number of significant digits in the output.
#' @param ... Additional arguments passed to the \code{wilcoxsign_test} 
#'            function.             
#' 
#' @details  r  is calculated as Z divided by 
#'           square root of the number of observations in one group. This 
#'           results in a statistic that ranges from -1 to 1.
#'           This range doesn't hold if \code{cases=FALSE}.
#'           
#'           This statistic reports a smaller effect size than does
#'           the matched-pairs rank biserial correlation coefficient 
#'           (\code{wilcoxonPairedRC}), and won't reach a value
#'           of -1 or 1 unless there are ties in paired differences.
#'
#'           Currently, the function makes no provisions for \code{NA}
#'           values in the data.  It is recommended that \code{NA}s be removed
#'           beforehand.
#'
#'           When the data in the first group are greater than
#'           in the second group, r is positive.
#'           When the data in the second group are greater than
#'           in the first group, r is negative.
#'           Be cautious with this interpretation, as R will alphabetize
#'           groups if \code{g} is not already a factor.
#'           
#'           When r is close to extremes,
#'           or with small counts in some cells,
#'           the confidence intervals 
#'           determined by this
#'           method may not be reliable, or the procedure may fail.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_06.html}
#' @seealso \code{\link{wilcoxonPairedRC}}
#' @concept effect size ordinal nominal
#' @return A single statistic, r.
#'         Or a small data frame consisting of r,
#'         and the lower and upper confidence limits.  
#'         
#' @examples
#' data(Pooh)
#' wilcox.test(Likert ~ Time, data=Pooh, paired=TRUE, exact=FALSE)
#' wilcoxonPairedR(x = Pooh$Likert, g = Pooh$Time)
#' 
#' @importFrom coin wilcoxsign_test
#' @importFrom boot boot boot.ci
#' 
#' @export
 
wilcoxonPairedR = function (x, g=NULL, coin=FALSE,
                            ci=FALSE, conf=0.95, type="perc",
                            R=1000, 
                            histogram=FALSE,
                            cases=TRUE,
                            digits=3, ... ){

  if(is.factor(g)==F){g=factor(g)}
  x = x[(as.numeric(g))<3]
  g = g[(as.numeric(g))<3]
  g = droplevels(g)
  
  if(coin){
    WT = suppressWarnings(wilcoxsign_test(x[as.numeric(g)==1] ~ x[as.numeric(g)==2], ...))
    Z  = as.numeric(statistic(WT, type="standardized"))
    }
  if(coin==FALSE){
    Z = wilcoxonZ(x = x[as.numeric(g)==1], y = x[as.numeric(g)==2], paired=TRUE)
    }
  n  = length(x[as.numeric(g)==1])
  if(cases==TRUE){r  = Z/sqrt(n)}
  if(cases==FALSE){r  = Z/sqrt(n*2)}
  RR = signif(r, digits=digits)
  
if(ci==TRUE){
  Data = data.frame(x1=x[as.numeric(g)==1], x2=x[as.numeric(g)==2])
  Function = function(input, index){
    Input = input[index,]
      if(coin){
        WT = suppressWarnings(wilcoxsign_test(x1 ~ x2,
                             data=Input, ...))
        Z  = as.numeric(statistic(WT, type="standardized"))
      }
      if(coin==FALSE){
        Z = wilcoxonZ(x = Input$x1, y = Input$x2, paired=TRUE)
      }
    n = length(Input$x1)
    if(cases==TRUE){r  = Z/sqrt(n)}
    if(cases==FALSE){r  = Z/sqrt(n*2)}
    RR = signif(r, digits=digits)
    return(RR)
    }
  Boot = boot(Data, Function, R=R)
  BCI  = boot.ci(Boot, conf=conf, type=type)
  if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3];}
  if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5];}
  if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5];}
  if(type=="bca") {CI1=BCI$bca[4];      CI2=BCI$bca[5];}  
  
  CI1=signif(CI1, digits=digits)
  CI2=signif(CI2, digits=digits)
  
  if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray", xlab="r", main="")}
}
  
if(ci==FALSE){names(RR)="r"; return(RR)}
if(ci==TRUE){DF=data.frame(r=RR, lower.ci=CI1, upper.ci=CI2)
             rownames(DF) = 1:nrow(DF)
             return(DF)
             }  
}  
