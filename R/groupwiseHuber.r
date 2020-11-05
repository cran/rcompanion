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
#' @param ... Other arguments passed to the \code{HuberM} function.
#'                
#' @details A wrapper for the \code{DescTools::HuberM} function
#'          to allow easy output for multiple groups.
#'          
#'          The input should include either \code{formula} and \code{data};
#'              or \code{data}, \code{var}, and \code{group}. (See examples). 
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
#' @references \url{http://rcompanion.org/rcompanion/d_08a.html}
#' @seealso \code{\link{groupwiseMean}}, \code{\link{groupwiseMedian}}
#'          \code{\link{groupwiseGeometric}}
#' @concept huber m-estimator confidence interval robust bootstrap 
#' @return A data frame of requested statistics by group.
#'          
#' @examples
#' ### Example with formula notation
#' data(Catbus)
#' groupwiseHuber(Steps ~ Teacher + Sex,
#'                data      = Catbus,
#'                ci.type   = "wald")
#'                
#' ### Example with variable notation
#' data(Catbus)
#' groupwiseHuber(data      = Catbus,
#'                var       = "Steps",
#'                group     = c("Teacher", "Sex"),
#'                ci.type   = "wald")
#'                                       
#' @importFrom DescTools HuberM
#' @importFrom plyr ddply rename
#' 
#' @export

groupwiseHuber = 
  function(formula=NULL, data=NULL,  var=NULL, group=NULL,
           conf.level=0.95, ci.type="wald", ...)
  {
 if(!is.null(formula)){
    var   = all.vars(formula[[2]])[1]
    group = all.vars(formula[[3]])
    }
  D1=
    ddply(.data=data,
          .variables=group, var,
          .fun=function(x, idx){
               sum(!is.na(x[,idx]))})
  funny = function(x, idx){HuberM(x[,idx], 
                           conf.level=conf.level, ci.type=ci.type, ...)}
  D2=
    ddply(.data=data,
          .variables=group, var,
          .fun=funny) 
                
D1 = rename(D1,c('V1'='n'))
D1$M.Huber = D2$hm
D1$lower.ci = D2$lwr.ci
D1$upper.ci = D2$upr.ci
return(D1)
}
