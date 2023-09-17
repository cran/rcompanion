#' @title Tukey's Ladder of Powers 
#'
#' @description Conducts Tukey's Ladder of Powers on a vector of values to 
#'              produce a more-normally distributed vector of values.
#' 
#' @param x A vector of values.
#' @param start The starting value of lambda to try.
#' @param end The ending value of lambda to try.
#' @param int The interval between lambda values to try.
#' @param statistic If \code{1}, uses Shapiro-Wilks test.
#'                  If \code{2}, uses Anderson-Darling test.
#' @param plotit If \code{TRUE}, produces plots of Shapiro-Wilks W or 
#'               Anderson-Darling A vs. lambda, a histogram of transformed
#'               values, and a quantile-quantile plot of transformed values.
#' @param verbose If \code{TRUE}, prints extra output for Shapiro-Wilks            
#'                W or Anderson-Darling A vs. lambda.
#' @param quiet If \code{TRUE}, doesn't print any output to the screen.
#' @param returnLambda If \code{TRUE}, returns only the lambda value,
#'                     not the vector of transformed values.
#'                
#' @details The function simply loops through lamdba values from \code{start}
#'          to \code{end} at an interval of \code{int}.
#'          
#'          The function then chooses the lambda which maximizes the 
#'          Shapiro-Wilks W statistic or minimizes the Anderson-Darling A
#'          statistic.
#'          
#'          It may be beneficial to add a constant to the input vector so that
#'          all values are posititive.  For left-skewed data, a (Constant - X)
#'          transformation may be helpful. Large values may need to be scaled. 
#'          
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{https://rcompanion.org/handbook/I_12.html}
#' 
#' @concept Tukey's ladder of powers
#' @concept normal distribution
#' @concept transformation
#' 
#' @return The transformed vector of values. The chosen lambda value is printed
#'         directly.
#'         
#' @examples
#' ### Log-normal distribution example
#' Conc = rlnorm(100)
#' Conc.trans = transformTukey(Conc)
#'                       
#' @importFrom stats rlnorm qqline qqnorm
#' @importFrom nortest ad.test
#' 
#' @export

transformTukey = 
   function(x, start=-10, end=10, int=0.025,
            plotit=TRUE, verbose=FALSE, quiet=FALSE, statistic=1, 
            returnLambda=FALSE)
   {
   n=(end-start)/int
   lambda=as.numeric(rep(0.00, n))
   W=as.numeric(rep(0.00, n))
   Shapiro.p.value=as.numeric(rep(0.00, n))
   if(statistic==2){
      A=as.numeric(rep(1000.00, n))
      Anderson.p.value=as.numeric(rep(0.00, n))
      }
   for(i in (1:n)){
      lambda[i] = signif(start+(i-1)*int, digits=4)
      if (lambda[i]>0) {TRANS = x ^ lambda[i]}
      if (lambda[i]==0){TRANS = log(x)}
      if (lambda[i]<0) {TRANS = -1 * x ^ lambda[i]}
      W[i]=NA
      if(statistic==2){A[i]=NA}
      if (any(is.infinite(TRANS))==FALSE & any(is.nan(TRANS))==FALSE)
         {
         W[i]=signif(shapiro.test(TRANS)$statistic, digits=4)
         Shapiro.p.value[i]=signif(shapiro.test(TRANS)$p.value, digits=4)
         if(statistic==2){
            A[i]=signif(ad.test(TRANS)$statistic, digits=4)
            Anderson.p.value[i]=signif(ad.test(TRANS)$p.value, digits=4)
         }
      } 
   }
   if(statistic==2){
      df=data.frame(lambda, W, Shapiro.p.value, A, Anderson.p.value)
   }
   if(statistic==1){
      df=data.frame(lambda, W, Shapiro.p.value)
   }
   if(verbose==TRUE){print(df)}
   if(plotit==TRUE){
      if(statistic==1){plot(lambda, W, col="black")}
      if(statistic==2){plot(lambda, A, col="blue")}
      }
   if(statistic==1){df2 = df[with(df, order(-W)),]}
   if(statistic==2){df2 = df[with(df, order(A)),]}
   if(quiet==FALSE){
   cat("\n")
   print(df2[1,])
   cat("\n")
   cat("if (lambda >  0){TRANS = x ^ lambda}","\n")
   cat("if (lambda == 0){TRANS = log(x)}","\n")
   cat("if (lambda <  0){TRANS = -1 * x ^ lambda}","\n")
   cat("\n")
   }
   lambda = df2[1,"lambda"]
   if (lambda>0) {TRANS = x ^ lambda}
   if (lambda==0){TRANS = log(x)}
   if (lambda<0) {TRANS = -1 * x ^ lambda}
   if(plotit==TRUE){
      plotNormalHistogram(TRANS, xlab="Transformed variable", 
                          linecol="red", col="lightgray")
      }
   if(plotit==TRUE){
   qqnorm(TRANS)
   qqline(TRANS, col="red")
   }
   if(returnLambda==FALSE){return(TRANS)}
   if(returnLambda==TRUE){
    names(lambda)="lambda"
    return(lambda)
   }
  }