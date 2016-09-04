#' @title Groupwise Huber M-estimators and confidence intervals 
#'
#' @description Calculates Huber M-estimator and confidence intervals for
#'              groups.
#' 
#' @param data The data frame to use.
#' @param var The measurement variable to use. The name is in double quotes.
#' @param group The grouping variable to use. The name is in double quotes.
#'              Multiple names are listed as a vector. (See example.)
#' @param conf.level The confidence interval to use.
#' @param conf.type The type of confidence interval to use. Can be
#'                  \code{"wald"} or \code{"boot"}.
#'                  See \code{\link{HuberM}} for details.
#' @param ... Other arguments passed to the \code{HuberM} function.
#'                
#' @details A wrapper for the \code{\link{HuberM}} function
#'          to allow easy output for multiple groups.
#'          
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/rcompanion/d_08a.html}
#' @seealso \code{\link{groupwiseMean}}, \code{\link{groupwiseMedian}}
#'          \code{\link{groupwiseGeometric}}
#' @concept huber m-estimator confidence interval robust bootstrap 
#' @return A data frame of requested statistics by group.
#'          
#' @examples
#' data(Catbus)
#' groupwiseHuber(data      = Catbus,
#'                var       = "Steps",
#'                group     = c("Teacher", "Sex"),
#'                conf.type = "wald")
#'                       
#' @importFrom DescTools HuberM
#' @importFrom plyr ddply rename
#' 
#' @export

groupwiseHuber = 
  function(data, group, var, conf.level=0.95, conf.type="wald", ...)
  {
  D1=
    ddply(.data=data,
          .variables=group, var,
          .fun=function(x, idx){
               sum(!is.na(x[,idx]))})
  funny = function(x, idx){HuberM(x[,idx], 
                           conf.level=conf.level, conf.type=conf.type, ...)}
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
