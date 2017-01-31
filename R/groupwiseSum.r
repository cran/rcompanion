#' @title Groupwise sums
#'
#' @description Calculates sums for groups.
#'
#' @param formula A formula indicating the measurement variable and
#'                the grouping variables. e.g. y ~ x1 + x2. 
#' @param data The data frame to use.
#' @param var The measurement variable to use. The name is in double quotes.
#' @param group The grouping variable to use. The name is in double quotes.
#'              Multiple names are listed as a vector. (See example.)
#' @param ... Other arguments passed to the \code{sum} function
#' @param digits The number of significant figures to use in output.
#' 
#' @details The input should include either \code{formula} and \code{data};
#'              or \code{data}, \code{var}, and \code{group}. (See examples).
#'                    
#' @note    The parsing of the formula is simplistic. The first variable on the
#'          left side is used as the measurement variable.  The variables on the
#'          right side are used for the grouping variables.
#'          
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @seealso \code{\link{groupwiseMean}}, \code{\link{groupwiseMedian}},
#'          \code{\link{groupwiseHuber}},
#'          \code{\link{groupwiseGeometric}}
#' @concept sum 
#' @return A data frame of statistics by group.
#'          
#' @examples
#' ### Example with formula notation
#' data(AndersonBias)
#' groupwiseSum(Count ~ Result + Sex,
#'              data        = AndersonBias)
#'                 
#' ### Example with variable notation
#' data(AndersonBias)
#' groupwiseSum(data        = AndersonBias,
#'              var         = "Count",
#'              group       = c("Result", "Sex"))
#'                       
#' @importFrom plyr ddply rename
#' 
#' @export

groupwiseSum = 
  function(formula=NULL, data=NULL, var=NULL, group=NULL,
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
  fun1 = function(x, idx){as.numeric(sum(x[,idx], na.rm=TRUE), ...)}
  D1=
     ddply(.data=data,
           .variables=group, var,
           .fun=fun1)
  ####################
DF = rename(DF,c('V1'='n'))
DF$Sum                                = signif(D1$V1, digits=digits)
return(DF)
}
