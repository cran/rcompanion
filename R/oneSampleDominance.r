#' @title Dominance statistic for one-sample data
#'
#' @description Calculates a dominance 
#'              effect size statistic
#'              compared with a theoretical median
#'              for one-sample data
#'              with confidence intervals by bootstrap
#'
#' @param x A vector of numeric values.
#' @param mu The median against which to compare the values.
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
#' @param na.rm If \code{TRUE}, removes \code{NA} values from
#'              the input vector \code{x}.
#' @param ... Additional arguments. 
#'             
#' @details The calculated \code{Dominance} statistic is simply 
#'          the proportion of observations greater than \code{mu} minus the
#'          the proportion of observations less than \code{mu}.
#'          
#'          It will range from -1 to 1, with 0 indicating that the median is
#'          equal to \code{mu}, 
#'          and 1 indicating that the observations are all greater in value
#'          than \code{mu},
#'          and -1 indicating that the observations are all less in value
#'          than \code{mu}.
#'          
#'          This statistic is appropriate for truly ordinal data,
#'          and could be considered an effect size statistic for
#'          a one-sample sign test. 
#'          
#'          Ordered category data need to re-coded as
#'          numeric, e.g. as with \code{as.numeric(Ordinal.variable)}.
#'
#'           When the statistic is close to 1 or close to -1,
#'           or with small sample size,
#'           the confidence intervals 
#'           determined by this
#'           method may not be reliable, or the procedure may fail.
#'           
#'           VDA is the analogous statistic, converted to a probability,
#'           ranging from 0 to 1, specifically,
#'           \code{VDA = Dominance / 2 + 0.5}.
#'                      
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{https://rcompanion.org/handbook/F_03.html}
#' 
#' @seealso \code{\link{pairedSampleDominance}},
#'          \code{\link{cliffDelta}},
#'          \code{\link{vda}}
#' 
#' @concept effect size
#' @concept dominance
#' @concept sign test
#' @concept confidence interval
#' 
#' @return A small data frame consisting of descriptive statistics,
#'         the dominance statistic,
#'         and potentially the lower and upper confidence limits.
#'                            
#' @examples
#' data(Catbus)
#' library(DescTools)
#' SignTest(Catbus$Rating, mu=5.5)
#' oneSampleDominance(Catbus$Rating, mu=5.5)
#' 
#' @importFrom boot boot boot.ci
#' 
#' @export

oneSampleDominance =
  function(x, mu=0,
           ci=FALSE, conf=0.95, type="perc", R=1000, histogram=FALSE,
           digits=3, na.rm=TRUE, ...)
  {
    
    if(na.rm){x = x[!is.na(x)]}
    
    Sum    = length(x)
    AdjSum = sum(x[!is.na(x)] != mu)
    
    if(AdjSum==0){
      Greater = 0
      Less    = 0
      Equal   = 1
    }
    
    if(AdjSum>0){
      Greater = sum(x>mu) / Sum
      Less    = sum(x<mu) / Sum
      Equal   = sum(x==mu) / Sum
    }
    
    es  = Greater - Less
    vda = es / 2 + 0.5
    
    if(ci==TRUE){
      y = rep(mu, Sum)
      Data = data.frame(x,y)
      Function = function(input, index){
        Input = input[index,]
        ess = (sum(Input$x>mu) - sum(Input$x<mu)) / length(Input$x)
        return(ess)
      }
      Boot = boot(Data, Function, R=R)
      BCI  = boot.ci(Boot, conf=conf, type=type)
      if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3];}
      if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5];}
      if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5];}
      if(type=="bca") {CI1=BCI$bca[4];      CI2=BCI$bca[5];}
      CI1=signif(CI1, digits=digits)
      CI2=signif(CI2, digits=digits)
      if(is.na(median(x))){CI1=NA; CI2=NA}
      CI3 = CI1 / 2 + 0.5
      CI4 = CI2 / 2 + 0.5
      
      if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray",
                               main="", xlab="dominance")}
      
      Out = data.frame(n=Sum, Median=median(x), mu=mu,
                       Less      = signif(Less, digits),
                       Equal     = signif(Equal, digits),
                       Greater   = signif(Greater, digits),
                       Dominance = signif(es, digits),
                       lower.ci  = CI1, 
                       upper.ci  = CI2,
                       VDA      =  signif(vda, digits),
                       lower.vda.ci  = CI3, 
                       upper.vda.ci  = CI4)
    }
    if(ci==FALSE){
      Out = data.frame(n=Sum, Median=median(x), mu=mu,
                       Less      = signif(Less, digits),
                       Equal     = signif(Equal, digits),
                       Greater   = signif(Greater, digits),
                       Dominance = signif(es, digits),
                       VDA      =  signif(vda, digits))
    }
    return(Out)
  }