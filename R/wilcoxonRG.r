#' @title Glass rank biserial correlation coefficient
#' 
#' @description Calculates Glass rank biserial correlation coefficient
#'              effect size
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
#' @param verbose If \code{TRUE}, prints information on factor levels and ranks.
#' @param na.last Passed to \code{rank}. For example, can be set to
#'                \code{TRUE} to assign \code{NA} values a minimum rank.
#' @param ... Additional arguments passed to \code{rank}             
#' 
#' @details  rg  is calculated as 2 times the difference of mean of ranks
#'           for each group divided by the total sample size.
#'           It appears that rg is equivalent to Cliff's delta.
#'  
#'           \code{NA} values can be handled by the \code{rank} function.
#'           In this case, using \code{verbose=TRUE} is helpful
#'           to understand how the \code{rg} statistic is calculated.
#'           Otherwise, it is recommended that \code{NA}s be removed
#'           beforehand.
#'           
#'           When the data in the first group are greater than
#'           in the second group, rg is positive.
#'           When the data in the second group are greater than
#'           in the first group, rg is negative.
#'           Be cautious with this interpretation, as R will alphabetize
#'           groups if \code{g} is not already a factor.
#'           
#'           When rg is close to extremes,
#'           or with small counts in some cells,
#'           the confidence intervals 
#'           determined by this
#'           method may not be reliable, or the procedure may fail.
#'                      
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references King, B.M., P.J. Rosopa, and E.W. Minium. 2011. 
#'             Statistical Reasoning in the Behavioral Sciences, 6th ed.
#' @seealso \code{\link{wilcoxonR}}
#' @concept correlation effect size ordinal nominal
#' @return A single statistic, rg.  
#'         Or a small data frame consisting of rg,
#'         and the lower and upper confidence limits.  
#'         
#' @examples
#' data(Breakfast)
#' Table = Breakfast[1:2,]
#' library(coin)
#' chisq_test(Table, scores = list("Breakfast" = c(-2, -1, 0, 1, 2)))
#' wilcoxonRG(Table)
#' 
#' data(Catbus)
#' wilcox.test(Steps ~ Sex, data = Catbus)
#' wilcoxonRG(x = Catbus$Steps, g = Catbus$Sex)
#' 
#' ### Example from King, Rosopa, and Minium
#' Criticism = c(-3, -2, 0, 0, 2, 5, 7, 9)
#' Praise = c(0, 2, 3, 4, 10, 12, 14, 19, 21)
#' Y = c(Criticism, Praise)
#' Group = factor(c(rep("Criticism", length(Criticism)),  
#'                 rep("Praise", length(Praise))))
#' wilcoxonRG(x = Y, g = Group, verbose=TRUE)
#' 
#' @importFrom boot boot boot.ci
#' 
#' @export
 
wilcoxonRG = function (x, g=NULL, group="row", 
                      ci=FALSE, conf=0.95, type="perc",
                      R=1000, histogram=FALSE, digits=3, 
                      reportIncomplete=FALSE,
                      verbose=FALSE, na.last=NA, ... ){
  
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
  
  RR=rank(x, na.last=na.last, ...)
  RRa = RR[as.numeric(g)==1]
  RRb = RR[as.numeric(g)==2]
  rG = 2 * (mean(RRa) - mean(RRb))/length(RR)
  rG=signif(rG, digits=digits)
  
  if(verbose)
  {
  cat("\n")
  cat("Levels:",levels(g))
  cat("\n")
  cat("n for",levels(g)[1], "=", length(RRa))
  cat("\n")
  cat("n for",levels(g)[2], "=",  length(RRb))
  cat("\n")
  cat("Mean of ranks for",levels(g)[1], "=",  mean(RRa))
  cat("\n")
  cat("Mean of ranks for",levels(g)[2], "=",  mean(RRb)) 
  cat("\n")
  cat("Difference in mean of ranks =",  (mean(RRa) - mean(RRb)))
  cat("\n")
  cat("Total n =",  length(RR))
  cat("\n")
  cat("2 * difference / total n =", rG)
  cat("\n")
  cat("\n")
  }
    
if(ci==TRUE){
  Data = data.frame(x,g)
  Function = function(input, index){
                    Input = input[index,]
                    if(length(unique(droplevels(Input$g)))==1){
                         FLAG=1
                         return(c(NA,FLAG))}
                    if(length(unique(droplevels(Input$g)))>1){
                         RR=rank(Input$x, na.last=na.last, ...)
                         RRa = RR[as.numeric(Input$g)==1]
                         RRb = RR[as.numeric(Input$g)==2]
                         RRa = RRa[!is.na(RRa)]
                         RRb = RRb[!is.na(RRb)]
                         rG = 2 * (mean(RRa) - mean(RRb))/length(RR)
                    FLAG=0  
                    return(c(rG, FLAG))}}
  Boot = boot(Data, Function, R=R)
  BCI  = boot.ci(Boot, conf=conf, type=type)
  if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3]}
  if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5]}
  if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5]}
  if(type=="bca")  {CI1=BCI$bca[4];     CI2=BCI$bca[5]}  
  
  if(sum(Boot$t[,2])>0 & reportIncomplete==FALSE) {CI1=NA; CI2=NA}
  
  CI1=signif(CI1, digits=digits)
  CI2=signif(CI2, digits=digits)
  
  if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray", xlab="rg",
                                       main="")}
  }
if(ci==FALSE){names(rG)="rg"; return(rG)}
if(ci==TRUE){DF=data.frame(rg=rG, lower.ci=CI1, upper.ci=CI2)
             rownames(DF) = 1:nrow(DF)
             return(DF)
             }  
}