#' @title r effect size for Wilcoxon two-sample rank-sum test
#' 
#' @description Calculates r effect size
#'              for Mann-Whitney, two-sample rank-sum test,
#'              or a table with an ordinal variable and a
#'              nominal variable with two levels; confidence intervals
#'              by bootstap.
#' 
#' @param x Either a two-way table or a two-way matrix.
#'          Can also be a vector of observations of an ordinal variable.
#' @param g If \code{x} is a vector, \code{g} is the vector of observations for
#'          the grouping, nominal variable.
#'          Only the first two levels of the nominal variable are used.
#' @param group If \code{x} is a table or matrix, \code{group} indicates whether
#'              the \code{"row"} or the \code{"column"} variable is
#'              the nominal, grouping variable.
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
#' @param ... Additional arguments passed to the \code{wilcox_test} function.             
#' 
#' @details  A Z value is extracted from the \code{wilcox_test} function in the
#'           coin package.  r  is calculated as Z divided by 
#'           square root of the total observations.
#'  
#'           Currently, the function makes no provisions for \code{NA}
#'           values in the data.  It is recommended that \code{NA}s be removed
#'           beforehand.
#'           
#'           When r is close to 0 or very large,
#'           or with small counts in some cells,
#'           the confidence intervals 
#'           determined by this
#'           method may not be reliable, or the procedure may fail.
#'                      
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_04.html}
#' @seealso \code{\link{freemanTheta}}
#' @concept correlation effect size ordinal nominal
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
#' data(PoohPiglet)
#' Data = PoohPiglet[PoohPiglet$Speaker %in% c("Pooh", "Piglet"),]
#' wilcox.test(Likert ~ Speaker, data = Data)
#' wilcoxonR(x = Data$Likert, g = Data$Speaker)
#' 
#' @importFrom coin wilcox_test
#' @importFrom boot boot boot.ci
#' 
#' @export
 
wilcoxonR = function (x, g=NULL, group="row", ci=FALSE, conf=0.95, type="perc",
                      R=1000, histogram=FALSE, digits=3, ... ){
  
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

  g = factor(g)
  x = x[as.numeric(g)<3]
  g = g[as.numeric(g)<3]
  g = droplevels(g)
  
  WT = suppressWarnings(wilcox_test(x ~ g, ...))
  Z  = as.numeric(statistic(WT, type="standardized"))
  N  = length(g)
  r  = abs(Z)/sqrt(N)
  RR = signif(r, digits=digits)

if(ci==TRUE){
  Data = data.frame(x,g)
  Function = function(input, index){
                    Input = input[index,]
                    WT = suppressWarnings(wilcox_test(x ~ g, data=Input, ...))
                    Z  = as.numeric(statistic(WT, type="standardized"))
                    N  = length(Input$g)
                    r  = abs(Z)/sqrt(N)
                    return(r)}
  Boot = boot(Data, Function, R=R)
  BCI  = boot.ci(Boot, conf=conf, type=type)
  if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3];}
  if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5];}
  if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5];}
  if(type=="bca") {CI1=BCI$bca[4];      CI2=BCI$bca[5];}  
  
  CI1=signif(CI1, digits=digits)
  CI2=signif(CI2, digits=digits)
  
  if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray")}

}
  
if(ci==FALSE){names(RR)="r"; return(RR)}
if(ci==TRUE){return(data.frame(r=RR, lower.ci=CI1, upper.ci=CI2))}  
  
}  
  
  