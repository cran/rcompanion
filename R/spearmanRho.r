#' @title Spearman's rho, Kendall's tau, Pearson's r
#'
#' @description Calculates Spearmans's rho, Kendall's tau, or Pearson's r
#'              with confidence intervals by bootstrap
#' 
#' @param formula A formula indicating the two paired variables,
#'                e.g.  \code{~ x + y}. The variables should be
#'                vectors of the same length.
#' @param data   The data frame to use. 
#' @param x If no formula is given, the values for one variable.
#' @param y The values for the other variable.
#' @param method One of "spearman", "kendall", or "pearson".
#'               Passed to \code{cor}.
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
#' @param ... Additional arguments passed to the \code{cor} function. 
#'             
#' @details This function is a wrapper for \code{stats::cor}
#'          with the addition of confidence intervals.
#'            
#'          The input should include either \code{formula} and \code{data};
#'          or \code{x}, and \code{y}.
#'          
#'           Currently, the function makes no provisions for \code{NA}
#'           values in the data.  It is recommended that \code{NA}s be removed
#'           beforehand.
#'
#'           When the returned statistic is close to -1 or close to 1,
#'           or with small sample size,
#'           the confidence intervals 
#'           determined by this
#'           method may not be reliable, or the procedure may fail.
#'                      
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/I_10.html}
#' @concept effect size
#' @return A single statistic, rho, tau, or r.
#'         Or a small data frame consisting of rho, tau, or r, 
#'         and the lower and upper confidence limits.
#'         
#' @examples
#' data(Catbus)
#' spearmanRho( ~ Steps + Rating, data=Catbus)
#' 
#' @importFrom stats cor
#' @importFrom boot boot boot.ci
#' 
#' @export

spearmanRho = 
 function(formula=NULL, data=NULL, x=NULL, y=NULL,
          method="spearman",
          ci=FALSE, conf=0.95, type="perc", R=1000, histogram=FALSE, digits=3,
          reportIncomplete=FALSE, ...){
   
  if(!is.null(formula)){
    x  = eval(parse(text=paste0("data","$",all.vars(formula)[1])))
    y  = eval(parse(text=paste0("data","$",all.vars(formula)[2])))
  }
  
  if(is.null(formula)){
   x = x
   y = y
  }
   
  Rho = signif(cor(x, y, method=method), digits=digits)
  
  if(ci==TRUE){
  Data = data.frame(x,y)
  Function = function(input, index){
                    Input = input[index,]
                    if(sd(Input$x)==0 | sd(Input$y)==0){
                       FLAG=1
                       return(c(NA,FLAG))
                       }
                    if(sd(Input$x)>0 & sd(Input$y)>0){
                      RHO = cor(Input$x, Input$y, method=method)
                      FLAG=0
                      return(c(RHO, FLAG))
                    }}
                    
  Boot = boot(Data, Function, R=R)
  
  BCI  = boot.ci(Boot, conf=conf, type=type)
  if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3];}
  if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5];}
  if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5];}
  if(type=="bca")  {CI1=BCI$bca[4];     CI2=BCI$bca[5];}
  
  if(sum(Boot$t[,2])>0 & reportIncomplete==FALSE) {CI1=NA; CI2=NA} 
  
  CI1=signif(CI1, digits=digits)
  CI2=signif(CI2, digits=digits)
  
  if(histogram==TRUE & method=="spearman"){
     hist(Boot$t[,1], col = "darkgray", main="", xlab="rho")}
  if(histogram==TRUE & method=="kendall"){
     hist(Boot$t[,1], col = "darkgray", main="", xlab="tau")}
  if(histogram==TRUE & method=="pearson"){
     hist(Boot$t[,1], col = "darkgray", main="", xlab="r")}
}
  
  if(ci==FALSE){names(Rho)="rho"}
  if(ci==FALSE & method=="kendall"){names(Rho)="tau"}
  if(ci==FALSE & method=="pearson"){names(Rho)="r"}
  if(ci==FALSE){return(Rho)}
  
  if(ci==TRUE){names(Rho) = ""}
  if(ci==TRUE & method=="spearman"){
     return(data.frame(rho=Rho, lower.ci=CI1, upper.ci=CI2))}
  if(ci==TRUE & method=="kendall"){
     return(data.frame(tau=Rho, lower.ci=CI1, upper.ci=CI2))}
  if(ci==TRUE & method=="pearson"){
     return(data.frame(r=Rho, lower.ci=CI1, upper.ci=CI2))}
 }