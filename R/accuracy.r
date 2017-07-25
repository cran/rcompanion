#' @title Minimum maximum accuracy, mean absolute percent error, and root mean square error
#'
#' @description Produces a table of fit statistics for multiple models.
#' 
#' @param fits A series of model object names. Must be a list.
#' @param plotit If \code{TRUE}, produces plots of the predicted values
#'               vs. the actual values for each model.
#' @param digits The number of significant digits in the output.              
#' @param ... Other arguments passed to \code{plot}.
#' 
#' @details  Produces a table of fit statistics for multiple models: 
#'           minimum maximum accuracy, mean absolute percentage error,
#'           root mean square error, normalized root mean square error,
#'           and accuracy based on normalized root mean square error.
#'           
#'           For minimum maximum accuracy, larger indicates
#'           a better fit, 
#'           and a perfect fit is equal to 1.
#'           
#'           For mean absolute percent error (MAPE), smaller
#'           indicates a better fit,
#'           and a perfect fit is equal to 0.
#'           
#'           Root mean square error (RMSE) has the same units as the predicted
#'           values.
#'           
#'           Normalized root mean square error (NRMSE) is RMSE divided by
#'           the mean or the median of the values of the dependent variable. 
#'           
#'           NRMSE accuracy values are calculated as 1 minus
#'           NRMSE.  Like minimum maximum accuracy, larger indicates a
#'           better fit, and a perfect fit is equal to 1.
#'           
#'           Model objects currently supported: lm, glm, nls, betareg, gls,
#'           lmer, lmerTest, rq, loess, gam, glm.nb, glmRob.
#'           
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/G_14.html}
#' @seealso \code{\link{compareLM}}, \code{\link{compareGLM}}, \code{\link{nagelkerke}}
#' @concept accuracy mape r-squared
#' @return A list of two objects: The series of model calls, and a data 
#'         frame of statistics for each model.
#'         
#' @examples
#' data(BrendonSmall)
#' BrendonSmall$Calories = as.numeric(BrendonSmall$Calories)
#' BrendonSmall$Calories2 = BrendonSmall$Calories ^ 2
#' model.1 = lm(Sodium ~ Calories, data = BrendonSmall)
#' model.2 = lm(Sodium ~ Calories + Calories2, data = BrendonSmall)
#' model.3 = glm(Sodium ~ Calories, data = BrendonSmall, family="Gamma")
#' quadplat = function(x, a, b, clx) {
#'           ifelse(x  < clx, a + b * x   + (-0.5*b/clx) * x   * x,
#'                            a + b * clx + (-0.5*b/clx) * clx * clx)}
#' model.4 = nls(Sodium ~ quadplat(Calories, a, b, clx),
#'               data = BrendonSmall,
#'               start = list(a=519, b=0.359, clx = 2300))
#' accuracy(list(model.1, model.2, model.3, model.4), plotit=FALSE)
#' 
#' @importFrom graphics plot
#' @importFrom stats predict residuals
#' 
#' @export

accuracy = 
function (fits, plotit=TRUE, digits=3, ...) 
{
 n = length(fits)
 Y = matrix(rep(NA,n),
            ncol=1)
 colnames(Y) = c("Call")
 rownames(Y) = seq(1,n)

 Z = data.frame(Min.max.accuracy=rep(NA,n),
                MAPE=rep(NA,n),
                MSE=rep(NA,n),
                RMSE=rep(NA,n),
                NRMSE.mean=rep(NA,n),
                NRMSE.median=rep(NA,n),
                NRMSE.mean.accuracy=rep(NA,n),
                NRMSE.median.accuracy=rep(NA,n),
                stringsAsFactors=FALSE) 
   for(i in 1:n)
    {
   TOGGLE=FALSE
   predy  = 1
   actual = 1
   if(class(fits[[i]])[1]=="lm"){predy  = predict(fits[[i]])
                                 actual = predict(fits[[i]]) +
                                          residuals(fits[[i]])
                                 call   = fits[[i]]$call
                                 TOGGLE = TRUE}
   if(class(fits[[i]])[1]=="glm"){predy  = predict(fits[[i]], type="response")
                                  actual = predict(fits[[i]], type="response") +
                                           residuals(fits[[i]], type="response")
                                  call   = fits[[i]]$call
                                  TOGGLE = TRUE}
   if(class(fits[[i]])[1]=="nls"){predy  = predict(fits[[i]])
                                  actual = predict(fits[[i]]) +
                                           residuals(fits[[i]])
                                  call   = fits[[i]]$call
                                  TOGGLE = TRUE} 
   if(class(fits[[i]])[1]=="betareg"){
                                  predy  = predict(fits[[i]], type="response")
                                  actual = predict(fits[[i]], type="response") +
                                           residuals(fits[[i]], type="response")
                                  call   = fits[[i]]$call
                                  TOGGLE = TRUE}
   if(class(fits[[i]])[1]=="gls"){
                                  predy  = predict(fits[[i]])
                                  actual = predict(fits[[i]]) +
                                           residuals(fits[[i]])
                                  call   = fits[[i]]$call
                                  TOGGLE = TRUE}
   if(class(fits[[i]])[1]=="lme"){
                                  predy  = predict(fits[[i]])
                                  actual = predict(fits[[i]]) +
                                           residuals(fits[[i]])
                                  call   = fits[[i]]$call
                                  TOGGLE = TRUE}
   if(class(fits[[i]])[1]=="lmerMod"){
                                  predy  = predict(fits[[i]])
                                  actual = predict(fits[[i]]) +
                                           residuals(fits[[i]])
                                  call   = fits[[i]]@call
                                  TOGGLE = TRUE}
   if(class(fits[[i]])[1]=="merModLmerTest"){
                                  predy  = predict(fits[[i]])
                                  actual = predict(fits[[i]]) +
                                           residuals(fits[[i]])
                                  call   = fits[[i]]@call
                                  TOGGLE = TRUE}
   if(class(fits[[i]])[1]=="rq"){
                                  predy  = predict(fits[[i]])
                                  actual = predict(fits[[i]]) +
                                           residuals(fits[[i]])
                                  call   = fits[[i]]$call
                                  TOGGLE = TRUE}
   if(class(fits[[i]])[1]=="loess"){
                                  predy  = predict(fits[[i]])
                                  actual = predict(fits[[i]]) +
                                           residuals(fits[[i]])
                                  call   = fits[[i]]$call
                                  TOGGLE = TRUE}
   if(class(fits[[i]])[1]=="gam"){
                                  predy  = predict(fits[[i]])
                                  actual = predict(fits[[i]]) +
                                           residuals(fits[[i]])
                                  call   = fits[[i]]$call
                                  TOGGLE = TRUE}
   if(class(fits[[i]])[1]=="negbin"){
                                  predy  = predict(fits[[i]], type="response")
                                  actual = predict(fits[[i]], type="response") +
                                           residuals(fits[[i]], type="response")
                                  call   = fits[[i]]$call
                                  TOGGLE = TRUE}
   if(class(fits[[i]])[1]=="glmRob"){
                                  predy  = predict(fits[[i]], type="response")
                                  actual = predict(fits[[i]], type="response") +
                                           residuals(fits[[i]], type="response")
                                  call   = fits[[i]]$call
                                  TOGGLE = TRUE}
   
    Y[i,] = "Not supported"

    if(TOGGLE==TRUE){Y[i,] = deparse(call)[1]}
      
     actual = unname(actual)
     predy  = unname(predy)
    
     data = data.frame(cbind(actual=actual, predy=predy))
     mma  = mean(apply(data, 1, min) / apply(data, 1, max))
     mape = mean(abs((predy - actual))/actual)
     mse  = mean((actual - predy)^2)
     rmse = sqrt(mse)
     nrmse_mean = rmse/mean(actual)
     nrmse_median = rmse/median(actual)
     acc_nrmse_mean = 1 - nrmse_mean
     acc_nrmse_median = 1 - nrmse_median
    
     Z[i,]=c(NA,NA)
     if(TOGGLE==TRUE){
     Z[i,]=c(signif(mma,              digits=digits),
             signif(mape,             digits=digits),
             signif(mse,              digits=digits),
             signif(rmse,             digits=digits),
             signif(nrmse_mean,       digits=digits),
             signif(nrmse_median,     digits=digits),
             signif(acc_nrmse_mean,   digits=digits),
             signif(acc_nrmse_median, digits=digits))
     if(plotit){plot(data$actual, data$predy, 
                      main=paste0("Model ", i),
                      xlab="Actual", ylab="Predicted",
                      ...)}
     }
     }     
   
 W = list(Y, Z)
 names(W) = c("Models",
              "Fit.criteria")
 return(W)            
 } 