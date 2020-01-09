#' @title Cliff's delta
#'
#' @description Calculates Cliff's delta
#'              with confidence intervals by bootstrap
#' 
#' @param formula A formula indicating the response variable and
#'                the independent variable. e.g. y ~ group.
#' @param data   The data frame to use. 
#' @param x If no formula is given, the response variable for one group.
#' @param y The response variable for the other group.
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
#' @param ... Additional arguments passed to the \code{wilcox.test} function. 
#'             
#' @details Cliff's delta is an effect size statistic appropriate
#'          in cases where a Wilcoxon-Mann-Whitney test might be used.
#'          It ranges from -1 to 1, with 0 indicating stochastic equality,
#'          and 1 indicating that the first group dominates the second.
#'          It is linearly related to Vargha and Delaney's A.
#'           
#'          The function calculates Cliff's delta from the "W" 
#'          U statistic from the
#'          \code{wilcox.test} function.
#'          Specifically, \code{VDA = U/(n1*n2); CD = (VDA-0.5)*2}.
#'            
#'          The input should include either \code{formula} and \code{data};
#'          or \code{x}, and \code{y}. If there are more than two groups,
#'          only the first two groups are used.
#'          
#'           Currently, the function makes no provisions for \code{NA}
#'           values in the data.  It is recommended that \code{NA}s be removed
#'           beforehand.
#'           
#'           When the data in the first group are greater than
#'           in the second group, Cliff's delta is positive.
#'           When the data in the second group are greater than
#'           in the first group, Cliff's delta is negative.
#'           Be cautious with this interpretation, as R will alphabetize
#'           groups in the formula interface if the grouping variable
#'           is not already a factor.
#'
#'           When Cliff's delta is close to 1 or close to -1,
#'           or with small sample size,
#'           the confidence intervals 
#'           determined by this
#'           method may not be reliable, or the procedure may fail.
#'                      
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_04.html}
#' @seealso \code{\link{vda}}, \code{\link{multiVDA}}
#' @concept effect size
#' @return A single statistic, Cliff's delta.
#'         Or a small data frame consisting of Cliff's delta,
#'         and the lower and upper confidence limits.
#'         
#' @note    The parsing of the formula is simplistic. 
#'          The first variable on the
#'          left side is used as the measurement variable.  
#'          The first variable on the
#'          right side is used for the grouping variable.   
#'          
#' @examples
#' data(Catbus)
#' cliffDelta(Steps ~ Sex, data=Catbus)
#' 
#' @importFrom stats wilcox.test
#' @importFrom boot boot boot.ci
#' 
#' @export

cliffDelta = 
 function(formula=NULL, data=NULL, x=NULL, y=NULL, 
          ci=FALSE, conf=0.95, type="perc", R=1000, histogram=FALSE, digits=3,
          reportIncomplete=FALSE,
          ...){

  if(!is.null(formula)){
    x  = eval(parse(text=paste0("data","$",all.vars(formula[[2]])[1])))
    g  = factor(eval(parse(text=paste0("data","$",all.vars(formula[[3]])[1]))))
    A  = x[g==levels(g)[1]]
    B  = x[g==levels(g)[2]]
  }
  
  if(is.null(formula)){
   A = x
   B = y
   x = c(A, B)
   g = factor(c(rep("A", length(A)), rep("B", length(B))))
  }
   
  n1  = length(A)
  n2  = length(B)
  U   = suppressWarnings(wilcox.test(x=A, y=B, ...))$statistic
  VDA = U / (n1 * n2) 
  CD = signif((VDA-0.5)*2, digits=digits)
  
  
  if(ci==TRUE){
  Data = data.frame(x,g)
  Function = function(input, index){
                    Input = input[index,]
                    if(length(unique(droplevels(Input$g)))==1){
                       FLAG=1
                       return(c(NA,FLAG))}  
                    if(length(unique(droplevels(Input$g)))>1){
                       U = suppressWarnings(wilcox.test(x ~ g, 
                                            data=Input, ...))$statistic
                       n1  = length(Input$x[Input$g==levels(Input$g)[1]])
                       n2  = length(Input$x[Input$g==levels(Input$g)[2]])
                       p   = U / (n1 * n2)
                       cd  = (p-0.5)*2
                       FLAG=0
                       return(c(cd, FLAG))}}
  Boot = boot(Data, Function, R=R)
  BCI  = boot.ci(Boot, conf=conf, type=type)
  if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3];}
  if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5];}
  if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5];}
  if(type=="bca") {CI1=BCI$bca[4];      CI2=BCI$bca[5];}
  
   if(sum(Boot$t[,2])>0 & reportIncomplete==FALSE) {CI1=NA; CI2=NA} 
  
  CI1=signif(CI1, digits=digits)
  CI2=signif(CI2, digits=digits)
  
  if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray",
                          main="", xlab="delta")}
}
  
  if(ci==FALSE){names(CD)="Cliff.delta"; return(CD)}
  if(ci==TRUE){names(CD) = ""
             return(data.frame(Cliff.delta=CD, lower.ci=CI1, upper.ci=CI2))}

 }