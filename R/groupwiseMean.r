#' @title Groupwise means and confidence intervals 
#'
#' @description Calculates means and confidence intervals for
#'              groups.
#' 
#' @param data The data frame to use.
#' @param var The measurement variable to use. The name is in double quotes.
#' @param group The grouping variable to use. The name is in double quotes.
#'              Multiple names are listed as a vector. (See example.)
#' @param conf The confidence interval to use.
#' @param R The number of bootstrap replicates to use for bootstrapped
#'          statistics.
#' @param boot If \code{TRUE}, includes the mean of the bootstrapped means.  
#'             This can be used as an estimate of the mean for
#'             the group.
#' @param traditional If \code{TRUE}, includes the traditional confidence
#'                    intervals for the group means, using the t-distribution.
#' @param normal If \code{TRUE}, includes the normal confidence
#'                    intervals for the group means by bootstrap.
#'                    See \code{\link{boot.ci}}.
#' @param basic If \code{TRUE}, includes the basic confidence
#'                    intervals for the group means by bootstrap.
#'                    See \code{\link{boot.ci}}.
#' @param percentile If \code{TRUE}, includes the percentile confidence
#'                    intervals for the group means by bootstrap.
#'                    See \code{\link{boot.ci}}.
#' @param bca If \code{TRUE}, includes the BCa confidence
#'                    intervals for the group means by bootstrap.
#'                    See \code{\link{boot.ci}}.                 
#' @param digits The number of significant figures to use in output.
#' @param ... Other arguments passed to the \code{boot} function.
#'                
#' @details With some options, the function may not handle missing values well.
#'          This seems to happen particularly with \code{bca = TRUE}.
#'          
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/C_03.html}
#' @seealso \code{\link{groupwiseMedian}}, \code{\link{groupwiseHuber}}
#'          \code{\link{groupwiseGeometric}}
#' @concept mean confidence interval bootstrap 
#' @return A data frame of requested statistics by group.
#'          
#' @examples
#' data(Catbus)
#' groupwiseMean(data        = Catbus,
#'               var         = "Steps",
#'               group       = c("Teacher", "Sex"),
#'               traditional = FALSE,
#'               percentile  = TRUE)
#'                       
#' @importFrom boot boot boot.ci
#' @importFrom plyr ddply rename
#' 
#' @export

groupwiseMean = 
  function(data, group, var, conf=0.95, R=5000,
           boot=FALSE, traditional=TRUE,
           normal=FALSE, basic=FALSE,
           percentile=FALSE, bca=FALSE, 
           digits=3, ...)
  {
  ####################
  DF=
    ddply(.data=data,
          .variables=group, var,
          .fun=function(x, idx){
               sum(!is.na(x[,idx]))})
  ####################
  fun1 = function(x, idx){as.numeric(mean(x[,idx], na.rm=TRUE))}
  D1=
     ddply(.data=data,
           .variables=group, var,
           .fun=fun1)
  ####################
  if(boot==TRUE){
  fun2 = function(x, idx){mean(boot(x[,idx],
                            function(y,j) mean(y[j]),
                            R=R, ...)$t[,1])}
  D2=ddply(.data=data,
           .variables=group, var,
           .fun=fun2)
  }
  ####################
  if(basic==TRUE){
  fun4 = function(x, idx){boot.ci(boot(x[,idx],
                            function(y,j) mean(y[j]),
                            R=R, ...), conf=conf, 
                                    type="basic", ...)$basic[4]}
  fun5 = function(x, idx){boot.ci(boot(x[,idx],
                            function(y,j) mean(y[j]),
                            R=R, ...), conf=conf, 
                                    type="basic", ...)$basic[5]}
  D4=ddply(.data=data,
         .variables=group, var,
         .fun=fun4)
  D5=ddply(.data=data,
         .variables=group, var,
         .fun=fun5)
  }
  ####################
  if(normal==TRUE){
     fun6 = function(x, idx){boot.ci(boot(x[,idx],
                                          function(y,j) mean(y[j]),
                                          R=R, ...), conf=conf, 
                                     type="norm", ...)$normal[2]}
     fun7 = function(x, idx){boot.ci(boot(x[,idx],
                                          function(y,j) mean(y[j]),
                                          R=R, ...), conf=conf, 
                                     type="norm", ...)$normal[3]}
     D6=ddply(.data=data,
              .variables=group, var,
              .fun=fun6)                    
     D7=ddply(.data=data,
              .variables=group, var,
              .fun=fun7)
  }
  ####################
  if(percentile==TRUE){
     fun8 = function(x, idx){boot.ci(boot(x[,idx],
                                          function(y,j) mean(y[j]),
                                          R=R, ...), conf=conf, 
                                     type="perc", ...)$percent[4]}
     fun9 = function(x, idx){boot.ci(boot(x[,idx],
                                          function(y,j) mean(y[j]),
                                          R=R, ...), conf=conf, 
                                     type="perc", ...)$percent[5]}
     D8=ddply(.data=data,
              .variables=group, var,
              .fun=fun8)                    
     D9=ddply(.data=data,
              .variables=group, var,
              .fun=fun9)
  }
  ####################
  if(bca==TRUE){
     fun10 = function(x, idx){boot.ci(boot(x[,idx],
                                          function(y,j) mean(y[j]),
                                          R=R, ...), conf=conf, 
                                     type="bca", ...)$bca[4]}
     fun11 = function(x, idx){boot.ci(boot(x[,idx],
                                          function(y,j) mean(y[j]),
                                          R=R, ...), conf=conf, 
                                     type="bca", ...)$bca[5]}
     D10=ddply(.data=data,
              .variables=group, var,
              .fun=fun10)                    
     D11=ddply(.data=data,
              .variables=group, var,
              .fun=fun11)
  }  
  ####################
  if(traditional==TRUE){
     Confy = function (x, ...){
        S = sd(x)
        N = length(x)
        Dist = conf + (1 - conf)/2
        Inty = qt(Dist, df = (N - 1)) * S/sqrt(N)
        return(Inty)
     }
     fun12 = function(x, idx){mean(x[,idx])-Confy(x[,idx])}
     fun13 = function(x, idx){mean(x[,idx])+Confy(x[,idx])}
     D12=ddply(.data=data,
               .variables=group, var,
               .fun=fun12)                    
     D13=ddply(.data=data,
               .variables=group, var,
               .fun=fun13)
  }  
  ####################  
DF = rename(DF,c('V1'='n'))
DF$Mean                                   = signif(D1$V1, digits=digits)
if(boot==TRUE){DF$Boot.mean               = signif(D2$V1, digits=digits)}
if(basic|normal|percentile|bca|traditional){DF$Conf.level = conf}
if(traditional==TRUE){DF$Trad.lower= signif(D12$V1, digits=digits)}
if(traditional==TRUE){DF$Trad.upper= signif(D13$V1, digits=digits)}
if(basic==TRUE){DF$Basic.lower            = signif(D4$V1, digits=digits)}
if(basic==TRUE){DF$Basic.upper            = signif(D5$V1, digits=digits)}
if(normal==TRUE){DF$Normal.lower          = signif(D6$V1, digits=digits)}
if(normal==TRUE){DF$Normal.upper          = signif(D7$V1, digits=digits)}
if(percentile==TRUE){DF$Percentile.lower  = signif(D8$V1, digits=digits)}
if(percentile==TRUE){DF$Percentile.upper  = signif(D9$V1, digits=digits)}
if(bca==TRUE){DF$Bca.lower                = signif(D10$V1, digits=digits)}
if(bca==TRUE){DF$Bca.upper                = signif(D11$V1, digits=digits)}

return(DF)
}
