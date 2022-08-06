#' @title Mangiafico's d
#'
#' @description Calculates Mangiafico's d, which is the difference in medians
#'              divided by the pooled median absolute deviation,
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
#' @param reportIncomplete If \code{FALSE} (the default),
#'                         \code{NA} will be reported in cases where there
#'                         are instances of the calculation of the statistic
#'                         failing during the bootstrap procedure.
#' @param verbose If \code{TRUE}, reports the median difference and MAD.
#' @param digits The number of significant digits in the output.      
#' @param ... Other arguments passed to \code{mad()}.
#'             
#' @details Mangiafico's d is an appropriate effect size statistic where
#'          Mood's median test, or another test comparing two medians,
#'          might be used.  Note that the response variable is treated
#'          as at least interval.
#'          
#'          For normal samples, the result will be somewhat similar to
#'          Cohen's d.
#'            
#'          The input should include either \code{formula} and \code{data};
#'          or \code{x}, and \code{y}. If there are more than two groups,
#'          only the first two groups are used.
#'          
#'          Currently, the function makes no provisions for \code{NA}
#'          values in the data.  It is recommended that \code{NA}s be removed
#'          beforehand.
#'           
#'          When the data in the first group are greater than
#'          in the second group, d is positive.
#'          When the data in the second group are greater than
#'          in the first group, d is negative.
#'          
#'          Be cautious with this interpretation, as R will alphabetize
#'          groups in the formula interface if the grouping variable
#'          is not already a factor.
#'
#'          When d is close to 0 or close to 1,
#'          or with small sample size,
#'          the confidence intervals 
#'          determined by this
#'          method may not be reliable, or the procedure may fail.
#'                      
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_05.html}
#' @seealso \code{\link{multiMangiaficoD}}
#' @concept effect size
#' @return A single statistic, d.
#'         Or a small data frame consisting of d,
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
#' mangiaficoD(Steps ~ Gender, data=Catbus)
#' 
#' @importFrom stats median mad
#' @importFrom boot boot boot.ci
#' 
#' @export

mangiaficoD = 
 function(formula=NULL, data=NULL, x=NULL, y=NULL, 
          ci=FALSE, conf=0.95, type="perc", R=1000, histogram=FALSE, 
          reportIncomplete=FALSE, verbose=FALSE, digits=3,
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
  
   MedianA = median(A)
   MedianB = median(B)
   MADA    = mad(A,...)
   MADB    = mad(B,...)
   DIFF    = MedianA - MedianB
   MAD     = sqrt((MADA^2 + MADB^2)/2)
   D       = DIFF / MAD
   
   if(verbose){
     G1        = levels(g)[1]
     G2        = levels(g)[2]
     Out = data.frame(
       Group     = c(G1, G2, "", G1, G2, ""),
       Statistic = c("Median", "Median", "Difference",
                     "MAD", "MAD", "Pooled MAD"),
       Value     = c(signif(MedianA, digits=digits),
                     signif(MedianB, digits=digits),
                     signif(DIFF, digits=digits),
                     signif(MADA, digits=digits),
                     signif(MADB, digits=digits),
                     signif(MAD, digits=digits))
     )
     cat("\n")
     print(Out)
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
                      A = Input$x[Input$g==levels(Input$g)[1]]
                      B = Input$x[Input$g==levels(Input$g)[2]]
                      MedianA = median(A)
                      MedianB = median(B)
                      MADA    = mad(A,...)
                      MADB    = mad(B,...)
                      DIFF    = MedianA - MedianB
                      MAD     = sqrt((MADA^2 + MADB^2)/2)
                      dd       = DIFF / MAD
                         FLAG=0}
                      
                      return(c(dd, FLAG))}
  
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
                      main="", xlab="Mangiafico's d")}
}
  
  if(ci==FALSE){names(D)="d"; return(signif(D, digits=digits))}
  if(ci==TRUE){names(D) = ""
             return(data.frame(d=signif(D, digits=digits),
                               lower.ci=CI1, upper.ci=CI2))}
 }