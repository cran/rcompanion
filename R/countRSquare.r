#' @title Count pseudo r-squared for logistic and other binary outcome models
#'
#' @description Produces the count pseudo 
#'              r-squared measure for models with a binary outcome.
#' 
#' @param fit The fitted model object for which to determine pseudo r-squared.
#'            \code{glm} and \code{glmmTMB} are supported.
#'            Others may work as well.
#' @param digits The number of digits in the outputted values.
#' @param suppressWarnings If \code{TRUE}, suppresses warning messages.
#' @param plotit If \code{TRUE}, produces a simple plot of
#'               actual vs. predicted values.
#' @param jitter If \code{TRUE}, jitters the "actual" values in the plot.
#' @param pch Passed to \code{plot}. 
#' @param ... Additional arguments.
#' 
#' @details  The count pseudo r-squared is simply the number of correctly
#'           predicted observations divided the total number of observations.
#'           
#'           This version is appropriate for models with a binary outcome.
#'           
#'           The adjusted value deducts the count of the most frequent
#'           outcome from both the numerator and the denominator.
#'           
#'           It is recommended that the model is fit on data in long
#'           format.  That is, that the \code{weight} option not be used in
#'           the model.
#'           
#'           The function makes no provisions for \code{NA} values.
#'           It is recommended that \code{NA} values be removed before
#'           the determination of the model.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{https://stats.oarc.ucla.edu/other/mult-pkg/faq/general/faq-what-are-pseudo-r-squareds/},
#'             \url{https://rcompanion.org/handbook/H_08.html},
#'             \url{https://rcompanion.org/rcompanion/e_06.html}
#'             
#' @seealso \code{\link{nagelkerke}}, 
#'          \code{\link{efronRSquared}}, 
#'          \code{\link{accuracy}}
#' 
#' @concept r squared
#' @concept pseudo r squared
#' @concept logistic regression
#' 
#' @return A list including a description of the submitted model,
#'         a data frame with the pseudo r-squared results,
#'         and a confusion matrix of the results.
#'          
#' @examples
#' data(AndersonBias)
#'
#' ### Covert data to long format
#'
#' Long = AndersonBias[rep(row.names(AndersonBias), AndersonBias$Count),
#'                     c("Result", "County", "Gender")]
#' rownames(Long) = seq(1:nrow(Long))
#' str(Long)
#'
#' ### Fit model and determine count r-square
#'
#' model = glm(Result ~ County + Gender + County:Gender,
#'             data = Long,
#'             family = binomial())
#'
#'countRSquare(model)
#' 
#' @importFrom stats predict addmargins
#' 
#' @export

countRSquare = 
function(fit, digits=3, suppressWarnings=TRUE, 
         plotit=FALSE, jitter=FALSE, pch=1, ...)
{ 
  if(suppressWarnings){
  Predy     = suppressWarnings(predict(fit, type="response"))
  Actual    = Predy + 
              suppressWarnings(residuals(fit, type = "response"))
  }
  if(!suppressWarnings){
    Predy     = predict(fit, type="response")
    Actual    = Predy + 
                residuals(fit, type = "response")
  }
  
  Correct   = sum(round(Predy)==Actual)
  Incorrect = sum(round(Predy)!=Actual)
  
  r2  = (Correct)/(Correct + Incorrect)
  k   = max(sum(Actual==0), sum(Actual==1))
  r2c = (Correct-k)/(Correct + Incorrect - k)
  
  if(plotit & jitter){
    plot(Predy, jitter(Actual), xlab="Predicted", ylab="Actual", pch=pch)
  }
  if(plotit & !jitter){
    plot(Predy, Actual, xlab="Predicted", ylab="Actual", pch=pch)
  }

  Y      = toString(fit$call)
  
  Z      = data.frame(Count.R2=NA, Count.R2.corrected=NA)
  Z[1,1] = signif(r2, digits=digits)
  Z[1,2] = signif(r2c, digits=digits)
  
  X = addmargins(table(Actual, round(Predy)))
  names(dimnames(X))=c("Actual", "Predicted")
  
  V = list(Y, Z, X) 
  names(V) = c("Model", "Result", "Confusion.matrix")
  return(V)            
}