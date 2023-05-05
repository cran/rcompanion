#' @title Grissom and Kim's Probability of Superiority (PS)
#'
#' @description Calculates Grissom and Kim's Probability of Superiority (PS)
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
#' @param verbose If \code{TRUE}, reports the proportion of ties and
#'                the proportions of (Ya > Yb) and (Ya < Yb).                         
#' @param ... Additional arguments, not used. 
#'             
#' @details PS is an effect size statistic appropriate
#'          in cases where a Wilcoxon-Mann-Whitney test might be used.
#'          It ranges from 0 to 1, with 0.5 indicating stochastic equality,
#'          and 1 indicating that the first group dominates the second.
#'          
#'          PS is defined as P(Ya > Yb), with no provision made for tied values
#'          across groups.
#'          
#'          If there are no tied values, PS will be equal to VDA. 
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
#'           in the second group, PS is greater than 0.5.
#'           When the data in the second group are greater than
#'           in the first group, PS is less than 0.5.
#'           
#'           Be cautious with this interpretation, as R will alphabetize
#'           groups in the formula interface if the grouping variable
#'           is not already a factor.
#'
#'           When PS is close to 0 or close to 1,
#'           or with small sample size,
#'           the confidence intervals 
#'           determined by this
#'           method may not be reliable, or the procedure may fail.
#'                      
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references Grissom, R.J. and J.J. Kim. 2012. Effect Sizes for Research. 
#'             2nd ed. Routledge, New York.
#'             
#'             \url{http://rcompanion.org/handbook/F_04.html}
#'             
#' @seealso \code{\link{cliffDelta}}, 
#'          \code{\link{vda}}
#'          
#' @concept effect size
#' @concept Grissom and Kim's probability of superiority
#' @concept Wilcoxon-Mann-Whitney
#' @concept confidence interval
#' 
#' @return A single statistic, PS.
#'         Or a small data frame consisting of PS,
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
#' wilcoxonPS(Steps ~ Gender, data=Catbus, verbose=TRUE)
#' 
#' @importFrom boot boot boot.ci
#' 
#' @export

wilcoxonPS = 
 function(formula=NULL, data=NULL, x=NULL, y=NULL, 
          ci=FALSE, conf=0.95, type="perc", R=1000, histogram=FALSE, digits=3,
          reportIncomplete=FALSE, verbose=FALSE,
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

    Matrix = outer(A,B,FUN="-")
    Diff   = Matrix>0
    PS    = signif(mean(Diff), digits=digits)
    
    if(verbose){
      Out = data.frame(
            Statistic = c("Proportion Ya > Yb","Proportion Ya < Yb",
                          "Proportion ties"),
            Value     = c(signif(mean(Matrix>0), digits=3),
                          signif(mean(Matrix<0), digits=3),
                          signif(mean(Matrix==0), digits=3))
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
                      
                      Matrix = outer(Input$x[Input$g==levels(Input$g)[1]],
                                     Input$x[Input$g==levels(Input$g)[2]],
                                     FUN="-")
                         Diff   = Matrix>0
                         p      = mean(Diff)
                         FLAG=0
                         return(c(p, FLAG))}}
  
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
                      main="", xlab="PS")}
}
  
  if(ci==FALSE){names(PS)="PS"; return(PS)}
  if(ci==TRUE){names(PS) = ""
             return(data.frame(PS=PS, lower.ci=CI1, upper.ci=CI2))}
 }