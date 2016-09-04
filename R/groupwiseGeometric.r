#' @title Groupwise geometric means and confidence intervals 
#'
#' @description Calculates geometric means and confidence intervals for
#'              groups
#' 
#' @param data The data frame to use.
#' @param var The measurement variable to use. The name is in double quotes.
#' @param group The grouping variable to use. The name is in double quotes.
#'              Multiple names are listed as a vector. (See example.)
#' @param conf The confidence interval to use.
#' @param na.rm If \code{TRUE}, removes NA values in the measurement variable.
#' @param digits The number of significant figures to use in output.
#' @param ... Other arguments.  Not currently useful.
#'                
#' @details The function computes means, standard deviations, standard errors, 
#'          and confidence intervals on log-transformed values. Confidence
#'          intervals are calculated in the traditional
#'          manner with the t-distribution.  These statistics assume that
#'          the data are log-normally distributed. For data not meeting this
#'          assumption, medians and confidence intervals by bootstrap may be more 
#'          appropriate. 
#'          
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/C_03.html}
#' @seealso \code{\link{groupwiseMean}}, \code{\link{groupwiseMedian}}
#'          \code{\link{groupwiseHuber}}
#' @concept geometric mean confidence interval standard deviation error 
#' @return A data frame of geometric means, standard deviations, standard 
#'         errors, and confidence intervals.
#'          
#' @examples
#' data(Catbus)
#' groupwiseGeometric(data   = Catbus,
#'                    var    = "Steps",
#'                    group  = c("Sex", "Teacher"))
#'                       
#' @importFrom stats sd qt na.omit
#' @importFrom plyr ddply rename
#' 
#' @export

groupwiseGeometric = function (data, group, var, conf=0.95, na.rm = TRUE, 
                                digits=3, ...) {
 Confy = function (x, ...){
        S = sd(x)
        N = length(x)
     Dist = conf + (1 - conf)/2
     Inty = qt(Dist, df = (N - 1)) * S/sqrt(N)
     return(Inty)
     }
 Geom = function (x, na.rm = na.rm) {
        if (na.rm) {x = na.omit(x)}
        return(exp(mean(log(x))))
        }
 SE = function (x, ...) {
      if (na.rm) {x = na.omit(x)}
      return(sd(x)/sqrt(length(x)))
      } 
 Gsd1 = function (x, na.rm = na.rm) {
        if (na.rm) {x <- na.omit(x)}
        return(exp(mean(log(x))-sd(log(x))))
        }
 Gsd2 = function (x, na.rm = na.rm) {
        if (na.rm) {x <- na.omit(x)}
        return(exp(mean(log(x))+sd(log(x))))
 }

 Gse1 = function (x, na.rm = na.rm) {
        if (na.rm) {x <- na.omit(x)}
       return(exp(mean(log(x))-SE(log(x))))
        }
 Gse2 = function (x, na.rm = na.rm) {
        if (na.rm) {x <- na.omit(x)}
        return(exp(mean(log(x))+SE(log(x))))
        }
 Gci1 = function (x, na.rm = na.rm) {
        if (na.rm) {x <- na.omit(x)}
        return(exp(mean(log(x)-Confy(log(x)))))
        }
 Gci2 = function (x, na.rm = na.rm) {
        if (na.rm) {x <- na.omit(x)}
        return(exp(mean(log(x)+Confy(log(x)))))
        }
      D1=
         ddply(.data=data,
               .variables=group, var,
               .fun=function(x, idx){
                  sum(!is.na(x[,idx]))})                  
                  
      fun2 = function(x, idx){Geom(x[,idx],na.rm=na.rm)}
      D2=
         ddply(.data=data,
               .variables=group, var,
               .fun=fun2)
      
      fun3 = function(x, idx){Gsd1(x[,idx],na.rm=na.rm)}
      D3=
         ddply(.data=data,
               .variables=group, var,
               .fun=fun3)
      
      fun4 = function(x, idx){Gsd2(x[,idx],na.rm=na.rm)}
      D4=
         ddply(.data=data,
               .variables=group, var,
               .fun=fun4)
      
      fun5 = function(x, idx){Gse1(x[,idx],na.rm=na.rm)}
      D5=
         ddply(.data=data,
               .variables=group, var,
               .fun=fun5)
      
      fun6 = function(x, idx){Gse2(x[,idx],na.rm=na.rm)}
      D6=
         ddply(.data=data,
               .variables=group, var,
               .fun=fun6)
      fun7 = function(x, idx){Gci1(x[,idx],na.rm=na.rm)}
      D7=
         ddply(.data=data,
               .variables=group, var,
               .fun=fun7)
      fun8 = function(x, idx){Gci2(x[,idx],na.rm=na.rm)}
      D8=
         ddply(.data=data,
               .variables=group, var,
               .fun=fun8)          
      
      D1 = rename(D1,c('V1'='n'))
      D1$Geo.mean = signif(D2$V1, digits=digits)
      D1$sd.lower = signif(D3$V1, digits=digits)
      D1$sd.upper = signif(D4$V1, digits=digits)
      D1$se.lower = signif(D5$V1, digits=digits)
      D1$se.upper = signif(D6$V1, digits=digits)
      D1$ci.lower = signif(D7$V1, digits=digits)
      D1$ci.upper = signif(D8$V1, digits=digits)
      
      return(D1)
   }