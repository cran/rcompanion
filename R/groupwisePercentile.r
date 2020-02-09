#' @title Groupwise percentiles and confidence intervals 
#'
#' @description Calculates percentiles and confidence intervals for
#'              groups.
#'
#' @param formula A formula indicating the measurement variable and
#'                the grouping variables. e.g. y ~ x1 + x2. 
#' @param data The data frame to use.
#' @param var If no formula is given, the measurement variable to use. 
#'            The name is in double quotes.
#' @param group The grouping variable to use. The name is in double quotes.
#'              Multiple names are listed as a vector. (See example.)
#' @param tau  The percentile to use, expressed as a quantile,
#'             e.g. 0.5 for median, 0.25 for 25th percentile.
#' @param type The \code{type} value passed to the \code{quantile} function.
#' @param conf The confidence interval to use.
#' @param R The number of bootstrap replicates to use for bootstrapped
#'          statistics.
#' @param boot If \code{TRUE}, includes the mean of the bootstrapped percentile.  
#'             This can be used as an estimate of the percentile for
#'             the group.
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
#' @details The input should include either \code{formula} and \code{data};
#'              or \code{data}, \code{var}, and \code{group}. (See examples).
#'           
#'          With some options, the function may not handle missing values well.
#'          This seems to happen particularly with \code{bca = TRUE}.
#'          
#' @note    The parsing of the formula is simplistic. The first variable on the
#'          left side is used as the measurement variable.  The variables on the
#'          right side are used for the grouping variables.
#'          
#'        Results for ungrouped (one-sample) data can be obtained by either
#'          setting the right side of the formula to 1, e.g.  y ~ 1, or by
#'          setting \code{group=NULL}.                
#'                    
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_15.html}
#' @seealso \code{\link{groupwiseMean}}, \code{\link{groupwiseHuber}}, 
#'          \code{\link{groupwiseGeometric}}, \code{\link{groupwiseMedian}}
#' @concept percentile confidence interval bootstrap 
#' @return A data frame of requested statistics by group
#'          
#' @examples
#' ### Example with formula notation
#' data(Catbus)
#' groupwisePercentile(Steps ~ Teacher + Sex,
#'                     data        = Catbus,
#'                     tau         = 0.25,
#'                     bca         = FALSE,
#'                     percentile  = TRUE,
#'                     R           = 1000)
#'                 
#' ### Example with variable notation
#' data(Catbus)
#' groupwisePercentile(data         = Catbus,
#'                     var         = "Steps",
#'                     group       = c("Teacher", "Sex"),
#'                     tau         = 0.25,
#'                     bca         = FALSE,
#'                     percentile  = TRUE,
#'                     R           = 1000)
#'                       
#' @importFrom boot boot boot.ci
#' @importFrom plyr ddply rename
#' @importFrom stats quantile
#' 
#' @export

groupwisePercentile = 
  function(formula=NULL, data=NULL, var=NULL, group=NULL,
           conf=0.95, tau=0.5, type=7, 
           R=5000,
           boot=FALSE,
           basic=FALSE, normal=FALSE,
           percentile=FALSE, bca=TRUE, 
           digits=3, ...)
  {
  if(!is.null(formula)){
    var   = all.vars(formula[[2]])[1]
    group = all.vars(formula[[3]])
    }
  ####################
  DF=
    ddply(.data=data,
          .variables=group, var,
          .fun=function(x, idx){
               sum(!is.na(x[,idx]))})
  ####################
  fun1 = function(x, idx){as.numeric(quantile(x[,idx], probs=tau, type=type, 
                                              na.rm=TRUE))}
  D1=
     ddply(.data=data,
           .variables=group, var,
           .fun=fun1)
  ####################
  if(boot==TRUE){
  fun2 = function(x, idx){mean(boot(x[,idx],
                            function(y,j) quantile(y[j], probs=tau, type=type),
                            R=R, ...)$t[,1])}
  D2=ddply(.data=data,
           .variables=group, var,
           .fun=fun2)
  }
  ####################
  if(basic==TRUE){
  fun4 = function(x, idx){boot.ci(boot(x[,idx],
                            function(y,j) quantile(y[j], probs=tau, type=type),
                            R=R, ...), conf=conf, 
                                    type="basic", ...)$basic[4]}
  fun5 = function(x, idx){boot.ci(boot(x[,idx],
                            function(y,j) quantile(y[j], probs=tau, type=type),
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
                                          function(y,j) quantile(y[j], probs=tau, type=type),
                                          R=R, ...), conf=conf, 
                                     type="norm", ...)$normal[2]}
     fun7 = function(x, idx){boot.ci(boot(x[,idx],
                                          function(y,j) quantile(y[j], probs=tau, type=type),
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
                                          function(y,j) quantile(y[j], probs=tau, type=type),
                                          R=R, ...), conf=conf, 
                                     type="perc", ...)$percent[4]}
     fun9 = function(x, idx){boot.ci(boot(x[,idx],
                                          function(y,j) quantile(y[j], probs=tau, type=type),
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
                                          function(y,j) quantile(y[j], probs=tau, type=type),
                                          R=R, ...), conf=conf, 
                                     type="bca", ...)$bca[4]}
     fun11 = function(x, idx){boot.ci(boot(x[,idx],
                                          function(y,j) quantile(y[j], probs=tau, type=type),
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
DF = rename(DF,c('V1'='n'))
DF$tau                                        = tau
DF$Percentile                                 = signif(D1$V1, digits=digits)
if(boot==TRUE){DF$Boot.percentile             = signif(D2$V1, digits=digits)}
if(basic|normal|percentile|bca){DF$Conf.level = conf}
if(basic==TRUE){DF$Basic.lower                = signif(D4$V1, digits=digits)}
if(basic==TRUE){DF$Basic.upper                = signif(D5$V1, digits=digits)}
if(normal==TRUE){DF$Normal.lower              = signif(D6$V1, digits=digits)}
if(normal==TRUE){DF$Normal.upper              = signif(D7$V1, digits=digits)}
if(percentile==TRUE){DF$Percentile.lower      = signif(D8$V1, digits=digits)}
if(percentile==TRUE){DF$Percentile.upper      = signif(D9$V1, digits=digits)}
if(bca==TRUE){DF$Bca.lower                    = signif(D10$V1, digits=digits)}
if(bca==TRUE){DF$Bca.upper                    = signif(D11$V1, digits=digits)}
return(DF)
}
