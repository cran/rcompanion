#' @title Epsilon-squared
#'
#' @description Calculates epsilon-squared for a table with one ordinal
#'              variable and one nominal variable; confidence intervals
#'              by bootstap.
#' 
#' @param x Either a two-way table or a two-way matrix.
#'          Can also be a vector of observations of an ordinal variable.
#' @param g If \code{x} is a vector, \code{g} is the vector of observations for
#'          the grouping, nominal variable.
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
#' @param ... Additional arguments passed to the \code{kruskal.test} function.             
#' 
#' @details  Epsilon-squared is used as a measure of association
#'           for the Kruskal-Wallis test or for a two-way
#'           table with one ordinal and one nominal variable.
#'
#'           Currently, the function makes no provisions for \code{NA}
#'           values in the data.  It is recommended that \code{NA}s be removed
#'           beforehand.
#'           
#'           When epsilon-squared is close to 0 or very large,
#'           or with small counts in some cells,
#'           the confidence intervals 
#'           determined by this
#'           method may not be reliable, or the procedure may fail.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/H_11.html}
#' @seealso \code{\link{freemanTheta}}
#' @concept correlation epsilon ordinal nominal
#' @return A single statistic, epsilon-squared.
#'         Or a small data frame consisting of epsilon-squared,
#'         and the lower and upper confidence limits.  
#'         
#' @examples
#' data(Breakfast)
#' library(coin)
#' chisq_test(Breakfast, scores = list("Breakfast" = c(-2, -1, 0, 1, 2)))
#' epsilonSquared(Breakfast)
#' 
#' data(PoohPiglet)
#' kruskal.test(Likert ~ Speaker, data = PoohPiglet)
#' epsilonSquared(x = PoohPiglet$Likert, g = PoohPiglet$Speaker)
#' 
#' ### Same data, as matrix of counts
#' data(PoohPiglet)
#' XT = xtabs( ~ Speaker + Likert , data = PoohPiglet)
#' epsilonSquared(XT)
#' 
#' @importFrom stats kruskal.test
#' @importFrom boot boot boot.ci
#' 
#' @export
 
epsilonSquared = function (x, g=NULL, group="row", ci=FALSE, conf=0.95, 
                           type="perc",
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

  g  = factor(g)
  g  = droplevels(g)
  n  = length(g)
  
  KW = kruskal.test(x, g, ...)

  e2 = KW$statistic / (n-1)
  E2 = signif(e2, digits=digits)
  
  
if(ci==TRUE){
  Data = data.frame(x,g)
  Function = function(input, index){
                      Input = input[index,]
                      n  = length(Input$g)
                      KW = kruskal.test(Input$x, Input$g)
                      e2 = KW$statistic / (n-1)
                    return(e2)}
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
  
if(ci==FALSE){names(E2) = "epsilon.squared"; return(E2)}
if(ci==TRUE){names(E2) = ""
             return(data.frame(epsilon.squared=E2, lower.ci=CI1, upper.ci=CI2))}  
}
