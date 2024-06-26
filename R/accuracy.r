#' @title Minimum maximum accuracy, mean absolute percent error,
#'        median absolute error, 
#'        root mean square error, coefficient of variation,
#'        and Efron's pseudo r-squared
#'
#' @description Produces a table of fit statistics for multiple models.
#' 
#' @param fits A series of model object names. 
#'             Must be a list of model objects or a single model object.
#' @param plotit If \code{TRUE}, produces plots of the predicted values
#'               vs. the actual values for each model.
#' @param digits The number of significant digits in the output.              
#' @param ... Other arguments passed to \code{plot}.
#' 
#' @details  Produces a table of fit statistics for multiple models: 
#'           minimum maximum accuracy, mean absolute percentage error,
#'           median absolute error,
#'           root mean square error, normalized root mean square error, 
#'           Efron's pseudo r-squared, and coefficient of variation.
#'           
#'           For minimum maximum accuracy, larger indicates
#'           a better fit, 
#'           and a perfect fit is equal to 1.
#'           
#'           For mean absolute error (MAE), smaller
#'           indicates a better fit,
#'           and a perfect fit is equal to 0.
#'           It has the same units as the dependent variable.
#'           Note that here, MAE is simply the mean of the absolute
#'           values of the differences of predicted values and the
#'           observed values 
#'           (\code{MAE = mean(abs(predy - actual))}).
#'           There are other definitions of MAE and similar-sounding
#'           terms.
#'           
#'           Median absolute error (MedAE) is similar, except employing
#'           the median rather than the mean.
#'           
#'           For mean absolute percent error (MAPE), smaller
#'           indicates a better fit,
#'           and a perfect fit is equal to 0. The result is reported
#'           as a fraction.  That is, a result of 0.1 is equal to 10 percent.
#'           
#'           Root mean square error (RMSE) has the same units as the predicted
#'           values.
#'           
#'           Normalized root mean square error (NRMSE) is RMSE divided by
#'           the mean or the median of the values of the dependent variable.
#'           
#'           Efron's pseudo r-squared is calculated as 1 minus the residual sum 
#'           of squares divided by the total sum of squares.  For linear models
#'           (\code{lm} model objects), Efron's pseudo r-squared will be equal  
#'           to r-squared.  For other models, it should not be interpreted
#'           as r-squared, but can still be useful as a relative measure.
#'           
#'           \code{CV.prcnt} is the coefficient of variation for the model.
#'           Here it is expressed as a percent.  That is, a result of 10 = 
#'           10 percent.
#'           
#'           Model objects currently supported: lm, glm, nls, betareg, gls,
#'           lme, lmer, lmerTest, glmmTMB, 
#'           rq, loess, gam, glm.nb, glmRob, mblm, and rlm.
#'           
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{https://rcompanion.org/handbook/G_14.html}
#' 
#' @seealso \code{\link{compareLM}}, 
#'          \code{\link{compareGLM}}, 
#'          \code{\link{nagelkerke}}
#'
#' @concept accuracy
#' @concept r squared
#' @concept pseudo r squared
#' @concept RMSE
#' @concept MAPE
#' @concept coefficient of variation
#' @concept CV
#' @concept MAE
#' 
#' @return A list of two objects: The series of model calls, and a data 
#'         frame of statistics for each model.
#'         
#' @examples
#' data(BrendonSmall)
#' BrendonSmall$Calories = as.numeric(BrendonSmall$Calories)
#' BrendonSmall$Calories2 = BrendonSmall$Calories ^ 2
#' model.1 = lm(Sodium ~ Calories, data = BrendonSmall)
#' 
#' accuracy(model.1, plotit=FALSE)
#' 
#' model.2 = lm(Sodium ~ Calories + Calories2, data = BrendonSmall)
#' model.3 = glm(Sodium ~ Calories, data = BrendonSmall, family="Gamma")
#' quadplat = function(x, a, b, clx) {
#'           ifelse(x  < clx, a + b * x   + (-0.5*b/clx) * x   * x,
#'                            a + b * clx + (-0.5*b/clx) * clx * clx)}
#' model.4 = nls(Sodium ~ quadplat(Calories, a, b, clx),
#'               data = BrendonSmall,
#'               start = list(a=519, b=0.359, clx = 2300))
#'               
#' accuracy(list(model.1, model.2, model.3, model.4), plotit=FALSE)
#' 
#' ### Perfect and poor model fits
#' X = c(1, 2,  3,  4,  5,  6, 7, 8, 9, 10, 11, 12)
#' Y = c(1, 2,  3,  4,  5,  6, 7, 8, 9, 10, 11, 12)
#' Z = c(1, 12, 13, 6, 10, 13, 4, 3, 5,  6, 10, 14)
#' perfect = lm(Y ~ X)
#' poor    = lm(Z ~ X)
#' accuracy(list(perfect, poor), plotit=FALSE)
#' 
#' @importFrom graphics plot
#' @importFrom stats predict residuals
#' 
#' @export

accuracy = 
function (fits, plotit=FALSE, digits=3, ...) 
{
 if(class(fits)[1]!="list"){fits=list(fits)}
  
 n = length(fits)
 Y = matrix(rep(NA,n),
            ncol=1)
 colnames(Y) = c("Call")
 rownames(Y) = seq(1,n)

 Z = data.frame(Min.max.accuracy=rep(NA,n),
                MAE=rep(NA,n),
                MedAE=rep(NA,n),
                MAPE=rep(NA,n),
                MSE=rep(NA,n),
                RMSE=rep(NA,n),
                NRMSE.mean=rep(NA,n),
                NRMSE.median=rep(NA,n),
                Efron.r.squared=rep(NA,n),
                CV.prcnt=rep(NA,n),
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
   if(class(fits[[i]])[1]=="lmerModLmerTest"){
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
   if(class(fits[[i]])[1]=="rlm"){
                                  predy  = predict(fits[[i]])
                                  actual = predict(fits[[i]]) + 
                                  residuals(fits[[i]])
                                  call   = fits[[i]]$call
                                  TOGGLE = TRUE}
   if(class(fits[[i]])[1]=="glmmTMB"){
                                  predy  = predict(fits[[i]], type="response")
                                  actual = predict(fits[[i]], type="response") + 
                                  residuals(fits[[i]], type="response")
                                  call   = fits[[i]]$call
                                  TOGGLE = TRUE}
   if(class(fits[[i]])[1]=="mblm"){
                                  predy  = predict(fits[[i]])
                                  actual = predict(fits[[i]]) + 
                                  residuals(fits[[i]])
                                  call   = fits[[i]]$call
                                  TOGGLE = TRUE}
   
    Y[i,] = "Not supported"

    if(TOGGLE==TRUE){Y[i,] = deparse(call)[1]}
      
     actual = unname(actual)
     predy  = unname(predy)
    
     data = data.frame(cbind(actual=actual, predy=predy))
     mma  = mean(apply(data, 1, min) / apply(data, 1, max))
     mae = mean(abs(predy - actual))
     medae = median(abs(predy - actual))
     mape = mean(abs((predy - actual)/actual))
     mse  = mean((actual - predy)^2)
     rmse = sqrt(mse)
     nrmse_mean = rmse/mean(actual)
     cv_prcnt      = nrmse_mean*100 
     nrmse_median = rmse/median(actual)
     Var = sum((actual - mean(actual))^2)
     RSS = sum((actual - predy)^2)
     var_r_squared = 1 - RSS / Var

     Z[i,]=rep(NA,10)
     if(TOGGLE==TRUE){
     Z[i,]=c(signif(mma,                digits=digits),
             signif(mae,                digits=digits),
             signif(medae,              digits=digits),
             signif(mape,               digits=digits),
             signif(mse,                digits=digits),
             signif(rmse,               digits=digits),
             signif(nrmse_mean,         digits=digits),
             signif(nrmse_median,       digits=digits),
             signif(var_r_squared,      digits=digits),
             signif(cv_prcnt,           digits=digits))
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