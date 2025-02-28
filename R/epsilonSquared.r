#' @title Epsilon-squared
#'
#' @description Calculates epsilon-squared 
#'              as an effect size statistic,
#'              following a Kruskal-Wallis test, 
#'              or for a table with one ordinal
#'              variable and one nominal variable; 
#'              confidence intervals by bootstrap
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
#' @param reportIncomplete If \code{FALSE} (the default),
#'                         \code{NA} will be reported in cases where there
#'                         are instances of the calculation of the statistic
#'                         failing during the bootstrap procedure.
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
#'           Because epsilon-squared is always positive, 
#'           if \code{type="perc"}, the confidence interval will
#'           never cross zero, and should not
#'           be used for statistical inference.
#'           However, if \code{type="norm"}, the confidence interval
#'           may cross zero. 
#'           
#'           When epsilon-squared is close to 0 or very large,
#'           or with small counts in some cells,
#'           the confidence intervals 
#'           determined by this
#'           method may not be reliable, or the procedure may fail.
#'           
#' @note     Note that epsilon-squared as calculated by this function
#'           is equivalent to the eta-squared, or r-squared, as
#'           determined by an anova on the rank-transformed values.
#'           Epsilon-squared for Kruskal-Wallis is typically
#'           defined this way in the literature.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references King, B.M., P.J. Rosopa, and E.W. Minium. 2018. 
#'             Statistical Reasoning in the Behavioral Sciences, 7th ed. Wiley.
#'             
#'             \url{https://rcompanion.org/handbook/F_08.html}
#'
#' @seealso \code{\link{multiVDA}},
#'          \code{\link{ordinalEtaSquared}}
#'             
#' @concept effect size
#' @concept Kruskal-Wallis
#' @concept epsilon squared
#' @concept confidence interval
#' 
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
                           R=1000, histogram=FALSE, digits=3, 
                           reportIncomplete=FALSE,
                           ... ){
  
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
                      if(length(unique(droplevels(Input$g)))==1){
                         FLAG=1
                         return(c(NA,FLAG))}  
                      if(length(unique(droplevels(Input$g)))>1){
                         KW = kruskal.test(Input$x, Input$g, ...)
                         e2 = KW$statistic / (n-1)
                         FLAG=0
                         return(c(e2, FLAG))
                    }}
  Boot = boot(Data, Function, R=R)
  
  BCI  = boot.ci(Boot, conf=conf, type=type)
  if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3]}
  if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5]}
  if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5]}
  if(type=="bca")  {CI1=BCI$bca[4];     CI2=BCI$bca[5]}  
  
  if(sum(Boot$t[,2])>0 & reportIncomplete==FALSE) {CI1=NA; CI2=NA} 
  
  CI1=signif(CI1, digits=digits)
  CI2=signif(CI2, digits=digits)
  
  if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray", xlab="epsilon-squared",
                                       main="")}
}
  
if(ci==FALSE){names(E2) = "epsilon.squared"; return(E2)}
if(ci==TRUE){names(E2) = ""
             return(data.frame(epsilon.squared=E2, lower.ci=CI1, upper.ci=CI2))}  
}
