#' @title Groupwise Huber M-estimators and confidence intervals 
#'
#' @description Calculates Huber M-estimator and confidence intervals for
#'              groups.
#'                      
#' @param formula A formula indicating the measurement variable and
#'                the grouping variables. e.g. y ~ x1 + x2.
#' @param data The data frame to use.
#' @param var The measurement variable to use. The name is in double quotes.
#' @param group The grouping variable to use. The name is in double quotes.
#'              Multiple names are listed as a vector. (See example.)
#' @param conf.level The confidence interval to use.
#' @param ci.type The type of confidence interval to use. Can be
#'                  \code{"wald"} or \code{"boot"}.
#'                  See \code{HuberM} for details.
#' @param digits The number of significant figures to use in output.
#' @param ... Other arguments passed to the \code{HuberM} function.
#'                
#' @details A wrapper for the \code{DescTools::HuberM} function
#'          to allow easy output for multiple groups.
#'          
#'          The input should include either \code{formula} and \code{data};
#'              or \code{data}, \code{var}, and \code{group}. (See examples).
#'              
#'          Results for ungrouped (one-sample) data can be obtained by either
#'          setting the right side of the formula to 1, e.g.  y ~ 1, or by
#'          setting \code{group=NULL}.  
#'          
#' @note    The parsing of the formula is simplistic. The first variable on the
#'          left side is used as the measurement variable.  The variables on the
#'          right side are used for the grouping variables.          
#'          
#'          It is recommended to remove \code{NA} values before using this
#'          function.  At the time of writing, \code{NA} values will cause the
#'          function to fail if confidence intervals are requested.
#'          
#'          At the time of writing, the \code{ci.type="boot"} option
#'          produces \code{NA} results. This is a result from the
#'          \code{DescTools::HuberM} function.            
#'              
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{http://rcompanion.org/rcompanion/d_08a.html}
#' 
#' @seealso \code{\link{groupwiseMean}}, 
#'          \code{\link{groupwiseMedian}},
#'          \code{\link{groupwiseGeometric}}
#'          
#' @concept summary statistics
#' @concept Huber M-estimator
#' @concept confidence interval
#' 
#' @return A data frame of requested statistics by group.
#'          
#' @examples
#' ### Example with formula notation
#' data(Catbus)
#' groupwiseHuber(Steps ~ Teacher + Gender,
#'                data      = Catbus,
#'                ci.type   = "wald")
#'                
#' ### Example with variable notation
#' data(Catbus)
#' groupwiseHuber(data      = Catbus,
#'                var       = "Steps",
#'                group     = c("Teacher", "Gender"),
#'                ci.type   = "wald")
#'
#' ### Example with NA value and without confidence intervals
#' data(Catbus)
#' Catbus1 = Catbus
#' Catbus1[1, 'Steps'] = NA
#' groupwiseHuber(Steps ~ Teacher + Gender,
#'                data      = Catbus1,
#'                conf.level   = NA)
#'                                       
#' @importFrom DescTools HuberM
#' @importFrom plyr ddply rename
#' 
#' @export

groupwiseHuber = 
  function(formula=NULL, data=NULL,  var=NULL, group=NULL,
           conf.level=0.95, ci.type="wald", digits=3, ...)
  {
 if(!is.null(formula)){
    var   = all.vars(formula[[2]])[1]
    group = all.vars(formula[[3]])
 }
    
  D1=
    ddply(.data=data,
          .variables=group, var,
          .fun=function(x, idx){
               length(x[,idx])})
  
  funny = function(x, idx){HuberM(x[,idx], 
                           conf.level=conf.level, ci.type=ci.type, ...)}
  
  D2=
    ddply(.data=data,
          .variables=group, var,
          .fun=funny)

D1 = rename(D1,c('V1'='n'))

if(is.na(conf.level)){
D1$M.Huber = signif(D2$V1, digits=digits)
}

if(!is.na(conf.level)){
D1$M.Huber = signif(D2$hm, digits=digits)
D1$lower.ci = signif(D2$lwr.ci, digits=digits)
D1$upper.ci = signif(D2$upr.ci, digits=digits)
}

return(D1)
}
