#' @title Groupwise quantiles 
#'
#' @description Calculates quantiles of an ordered factor or numeric variable
#'              by one or more grouping variables 
#'          
#' @param formula A formula indicating the measurement variable and
#'                the grouping variables. e.g. y ~ x1 + x2.
#' @param data The data frame to use.
#' @param var The measurement variable to use. The name is in double quotes.
#' @param group The grouping variable to use. The name is in double quotes.
#'              Multiple names are listed as a vector. (See example.)
#' @param type The type of quantiles, as in the \code{quantile} function.
#'             For ordered factors, must be 1 or 3.
#' @param quantiles A vector of quantiles to determine.
#' @param digits For numeric variables, 
#'               the number of significant figures to use in output.
#' @param ... Other arguments passed to the \code{quantile} function.
#' 
#' @details The input should include either \code{formula} and \code{data};
#'              or \code{data}, \code{var}, and \code{group}. (See examples).
#'          
#'          Results for ungrouped (one-sample) data can be obtained by either
#'          setting the right side of the formula to 1, e.g.  y ~ 1, or by
#'          setting \code{group=NULL} when using \code{var}.
#'          
#' @note    The parsing of the formula is simplistic. The first variable on the
#'          left side is used as the measurement variable.  The variables on the
#'          right side are used for the grouping variables.
#'          
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{https://rcompanion.org/handbook/C_02.html}
#' 
#' @seealso \code{\link{groupwiseMedian}}, 
#'          \code{\link{groupwisePercentile}}
#'          
#' @concept summary statistics
#' @concept quantile
#' @concept percentile
#' 
#' @return A data frame of requested statistics by group.
#'          
#' @examples
#' Education = c("BA","PHD","BA","MA","HS","MA","HS","BA","BA","BA","MA","MA",
#'               "AA","AA","BA","BA","PHD")
#' 
#' Education = factor(Education, ordered=TRUE, levels=c("HS","AA","BA","MA","PHD"))
#' 
#' Gender = c("female","female","male","female","male","female","male","female",
#'            "female","female","male","male","female","male","other","female",
#'            "female")
#' 
#' County = c("Elwood","Bear Lake","Elwood","Elwood","Elwood","Bear Lake","Elwood",
#'            "Elwood","Elwood","Elwood","Elwood","Bear Lake","Elwood","Elwood",
#'            "Elwood","Elwood","Bear Lake")
#' 
#' Arthur = data.frame(Gender, County, Education)
#' 
#' Arthur$Ed.code = as.numeric(Arthur$Education)
#' 
#' ### One-way, two-way, and one-sample for ordered factor data
#' 
#' groupwiseFiveNumber(Education ~ Gender, data=Arthur, type=3)
#' 
#' groupwiseFiveNumber(Education ~ Gender + County, data=Arthur, type=3)
#' 
#' groupwiseFiveNumber(Education ~ 1, data=Arthur, type=3)
#' 
#' ### For numeric data
#' 
#' groupwiseFiveNumber(Ed.code ~ Gender, data=Arthur)
#' 
#' ### Variable notation
#' 
#' groupwiseFiveNumber(data=Arthur, var="Education", group="Gender", type=3)
#' 
#' ### Asking for different quantiles
#' 
#' groupwiseFiveNumber(Education ~ Gender, data=Arthur, type=3,
#'                     quantiles=c(0.00, 0.05, 0.10, 0.15, 0.20, 0.25))
#' 
#' ### Handling missing data
#' 
#' Arthur2 = Arthur
#' 
#' Arthur2$Education[1] = NA
#' Arthur2$Education[2] = NA
#' 
#' groupwiseFiveNumber(Education ~ Gender, data=Arthur2, type=3)
#'                                      
#' @importFrom plyr ddply rename
#' @importFrom stats quantile
#' 
#' @export

groupwiseFiveNumber = 
  function(formula=NULL, data=NULL,  var=NULL, group=NULL,
           type=7,
           quantiles=c(0.00, 0.25, 0.50, 0.75, 1.00),
           digits=3, ...)
  {
 if(!is.null(formula)){
    var   = all.vars(formula[[2]])[1]
    group = all.vars(formula[[3]])
 }

 nValid = NULL
   
  D1=
    ddply(.data=data,
          .variables=group, var,
          .fun=function(x, idx){
               length(x[,idx])})
  
  funny = function(x, idx){quantile(x[,idx], 
                           probs = quantiles,
                           type=type, na.rm=TRUE, ...)}
  
  D1 = rename(D1,c('V1'='n'))
 
  D2=
   ddply(.data=data,
         .variables=group, var,
         .fun=function(x, idx){sum(!is.na(x[,idx]))})
  
  D2 = rename(D2,c('V1'='nValid'))
  
  D3=
    ddply(.data=data,
          .variables=group, var,
          .fun=funny)

D3 = D3[, (ncol(D3) - length(quantiles) + 1) : ncol(D3)] 

D1$nValid = D2$nValid

if(is.numeric(data[[var]]) | is.integer(data[[var]])){
 for(i in 1 : ncol(D3)){
  D3[,i] = signif(D3[,i], digits)
 }
}

D1 = if(sum(D1$n != D1$nValid) == 0) {subset(D1, select=-c(nValid))} 
     else {D1}

D1 = cbind(D1, D3)

return(D1)
}
