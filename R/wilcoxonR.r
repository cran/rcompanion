#' @title r effect size for Wilcoxon two-sample rank-sum test
#' 
#' @description Calculates r effect size
#'              for Mann-Whitney two-sample rank-sum test,
#'              or a table with an ordinal variable and a
#'              nominal variable with two levels; confidence intervals
#'              by bootstrap.
#' 
#' @param x Either a two-way table or a two-way matrix.
#'          Can also be a vector of observations.
#' @param g If \code{x} is a vector, \code{g} is the vector of observations for
#'          the grouping, nominal variable.
#'          Only the first two levels of the nominal variable are used.
#' @param group If \code{x} is a table or matrix, \code{group} indicates whether
#'              the \code{"row"} or the \code{"column"} variable is
#'              the nominal, grouping variable.
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
#' @param reportIncomplete If \code{FALSE} (the default),
#'                         \code{NA} will be reported in cases where there
#'                         are instances of the calculation of the statistic
#'                         failing during the bootstrap procedure.
#' @param digits The number of significant digits in the output.
#' @param ... Additional arguments passed to the \code{wilcox_test} function.             
#' 
#' @details  r  is calculated as Z divided by 
#'           square root of the total observations.
#'           
#'           This statistic reports a smaller effect size than does
#'           Glass rank biserial correlation coefficient 
#'           (\code{wilcoxonRG}), and cannot reach 
#'           -1 or 1.  This effect is exaserbated when sample sizes
#'           are not equal.
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
#' 
#' @references \url{https://rcompanion.org/handbook/F_04.html}
#' 
#' @seealso \code{\link{freemanTheta}}, 
#'          \code{\link{wilcoxonRG}}
#'          
#' @concept effect size
#' @concept Wilcoxon-Mann-Whitney
#' @concept confidence interval
#' 
#' @return A single statistic, r.  
#'         Or a small data frame consisting of r,
#'         and the lower and upper confidence limits.  
#'         
#' @examples
#' data(Breakfast)
#' Table = Breakfast[1:2,]
#' library(coin)
#' chisq_test(Table, scores = list("Breakfast" = c(-2, -1, 0, 1, 2)))
#' wilcoxonR(Table)
#' 
#' data(Catbus)
#' wilcox.test(Steps ~ Gender, data = Catbus)
#' wilcoxonR(x = Catbus$Steps, g = Catbus$Gender)
#' 
#' @importFrom coin wilcox_test
#' @importFrom boot boot boot.ci
#' 
#' @export
 
wilcoxonR = function (x, g=NULL, group="row", coin=FALSE, 
                      ci=FALSE, conf=0.95, type="perc",
                      R=1000, histogram=FALSE, digits=3, 
                      reportIncomplete=FALSE, ... ){
  
  if(is.matrix(x)){x=as.table(x)}
  
  if(is.table(x)){
    Counts = as.data.frame(x)
    Long = Counts[rep(row.names(Counts), Counts$Freq), c(1, 2)]
    rownames(Long) = seq(1:nrow(Long))
    if(group=="row"){
       g=factor(Long[,1])
       x=as.numeric(Long[,2])}
    if(group=="column"){
       g=factor(Long[,2])
       x=as.numeric(Long[,1])}
  }

  if(is.factor(g)==F){g=factor(g)}
  x = x[as.numeric(g)<3]
  g = g[as.numeric(g)<3]
  g = droplevels(g)
  
  if(coin){
    WT = suppressWarnings(wilcox_test(x ~ g, ...))
    Z  = as.numeric(statistic(WT, type="standardized"))
  }
  if(coin==FALSE){
     Z = wilcoxonZ(x = x[as.numeric(g)==1], y = x[as.numeric(g)==2])
  }
  N  = length(g)
  r  = Z/sqrt(N)
  RR = signif(r, digits=digits)

if(ci==TRUE){
  Data = data.frame(x,g)
  Function = function(input, index){
                    Input = input[index,]
                    if(length(unique(droplevels(Input$g)))==1){
                         FLAG=1
                         return(c(NA,FLAG))}
                    if(length(unique(droplevels(Input$g)))>1){
                    if(coin){
                      WT = suppressWarnings(wilcox_test(x ~ g, data=Input, ...))
                      Z  = as.numeric(statistic(WT, type="standardized"))
                      }
                    if(coin==FALSE){
                      Z = wilcoxonZ(x = Input$x[as.numeric(Input$g)==1], 
                                    y = Input$x[as.numeric(Input$g)==2])
                      }
                    N  = length(Input$g)
                    r  = Z/sqrt(N)
                    FLAG=0
                    return(c(r, FLAG))}}
  Boot = boot(Data, Function, R=R)
  BCI  = boot.ci(Boot, conf=conf, type=type)
  if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3]}
  if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5]}
  if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5]}
  if(type=="bca")  {CI1=BCI$bca[4];     CI2=BCI$bca[5]}  
  
  if(sum(Boot$t[,2])>0 & reportIncomplete==FALSE) {CI1=NA; CI2=NA}
  
  CI1=signif(CI1, digits=digits)
  CI2=signif(CI2, digits=digits)
  
  if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray", xlab="r",
                                       main="")}
  }
if(ci==FALSE){names(RR)="r"; return(RR)}
if(ci==TRUE){DF=data.frame(r=RR, lower.ci=CI1, upper.ci=CI2)
             rownames(DF) = 1:nrow(DF)
             return(DF)
             }  
}
  
  