#' @title Efron's pseudo r-squared
#'
#' @description Produces Efron's pseudo r-squared from certain models,
#'              or vectors of
#'              residuals, predicted values, and actual values.
#'              Alternately produces minimum maximum accuracy, 
#'              mean absolute percent error, 
#'              root mean square error, or coefficient of variation.
#' 
#' @param model A model of the class lm, glm, nls, betareg, gls,
#'                                   lme, lmerMod, lmerModLmerTest, glmmTMB,
#'                                   rq, loess, gam, negbin, glmRob, rlm,
#'                                   or mblm.
#'              
#' @param actual A vector of actual y values
#' @param residual A vector of residuals
#' @param predicted A vector of predicted values
#' @param statistic The statistic to produce.
#'                  One of \code{"EfronRSquared"},
#'                         \code{"MinMaxAccuracy"},
#'                         \code{"MAE"},
#'                         \code{"MAPE"},
#'                         \code{"MSE"},
#'                         \code{"RMSE"},
#'                         \code{"NRMSE.Mean"},
#'                         \code{"CV"}.
#' @param plotit If \code{TRUE}, produces plots of the predicted values
#'               vs. the actual values.
#' @param digits The number of significant digits in the output.              
#' @param ... Other arguments passed to \code{plot}.
#' 
#' @details  Efron's pseudo r-squared is calculated as 1 minus the residual sum 
#'           of squares divided by the total sum of squares.  For linear models
#'           (\code{lm} model objects), Efron's pseudo r-squared will be equal  
#'           to r-squared.
#'           
#'           This function produces the same statistics as does the
#'           \code{accuracy} function.
#'           While the \code{accuracy} function extracts values from a model
#'           object, this function allows for the manual entry
#'           of residual, predicted, or actual values.
#'           
#'           It is recommended that the user consults the \code{accuracy}
#'           function
#'           for further details on these statistics, such as if the reported
#'           value is presented as a percentage or fraction.
#'           
#'           If \code{model}is not supplied,
#'           two of the following need to passed to the function:
#'           \code{actual}, \code{predicted}, \code{residual}.
#'
#'           Note that, for some model objects, to extract residuals
#'           and predicted values on the original scale,
#'           a \code{type="response"}
#'           option needs to be added to the call, e.g.
#'           \code{residuals(model.object, type="response")}.
#'
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{https://rcompanion.org/handbook/F_16.html}
#' 
#' @seealso \code{\link{accuracy}}, 
#'          \code{\link{nagelkerke}}
#'
#' @concept pseudo r squared
#' @concept r squared
#' @concept accuracy
#' @concept RMSE
#' @concept MAPE
#' @concept coefficient of variation
#' @concept MAE
#' 
#' @return A single statistic
#'         
#' @examples
#' data(BrendonSmall)
#' BrendonSmall$Calories = as.numeric(BrendonSmall$Calories)
#' BrendonSmall$Calories2 = BrendonSmall$Calories ^ 2
#' model.1 = lm(Sodium ~ Calories + Calories2, data = BrendonSmall)
#' 
#' efronRSquared(model.1)
#' 
#' efronRSquared(model.1, statistic="MAPE")
#' 
#' efronRSquared(actual=BrendonSmall$Sodium, residual=model.1$residuals)
#' efronRSquared(residual=model.1$residuals, predicted=model.1$fitted.values)
#' efronRSquared(actual=BrendonSmall$Sodium, predicted=model.1$fitted.values)
#' 
#' summary(model.1)$r.squared
#' 
#' @importFrom graphics plot
#' 
#' @export

efronRSquared = 
function (model=NULL, actual=NULL, predicted=NULL, residual=NULL, 
          statistic="EfronRSquared",
          plotit=FALSE, digits=3, ...) 
{
 if(!is.null(model)){
   predicted  = as.numeric(NA)
   residual   = as.numeric(NA)
    if(class(model)[1]=="lm"){
                              predicted  = predict(model)
                              residual = residuals(model)}
   if(class(model)[1]=="glm"){
                              predicted  = predict(model, type="response")
                              residual = residuals(model, type="response")}
   if(class(model)[1]=="nls"){
                              predicted  = predict(model)
                              residual = residuals(model)} 
   if(class(model)[1]=="betareg"){
                              predicted  = predict(model, type="response")
                              residual = residuals(model, type="response")}
   if(class(model)[1]=="gls"){predicted  = predict(model)
                              residual = residuals(model)}
   if(class(model)[1]=="lme"){
                              predicted  = predict(model)
                              residual = residuals(model)}
   if(class(model)[1]=="lmerMod"){
                              predicted  = predict(model)
                              residual = residuals(model)}
   if(class(model)[1]=="lmerModLmerTest"){
                              predicted  = predict(model)
                              residual = residuals(model)}
   if(class(model)[1]=="rq"){
                              predicted  = predict(model)
                              residual = residuals(model)}
   if(class(model)[1]=="loess"){
                              predicted  = predict(model)
                              residual = residuals(model)}
   if(class(model)[1]=="gam"){
                              predicted  = predict(model)
                              residual = residuals(model)}
   if(class(model)[1]=="negbin"){
                              predicted  = predict(model, type="response")
                              residual = residuals(model, type="response")}
   if(class(model)[1]=="glmRob"){
                              predicted  = predict(model, type="response")
                              residual = residuals(model, type="response")}
   if(class(model)[1]=="rlm"){
                              predicted  = predict(model)
                              residual = residuals(model)}
   if(class(model)[1]=="glmmTMB"){
                               predicted  = predict(model, type="response")
                               residual = residuals(model, type="response")}
   if(class(model)[1]=="mblm"){
     predicted  = predict(model)
     residual = residuals(model)}
      }

  if(is.null(actual)){    actual    = predicted + residual}
  if(is.null(residual)){  residual  = actual    - predicted}
  if(is.null(predicted)){ predicted = actual    - residual}
  
  actual     = unname(actual)
  residual   = unname(residual)
  predicted  = unname(predicted)
    
  data = data.frame(cbind(actual=actual, predicted=predicted))
  
  mma           = mean(apply(data, 1, min) / apply(data, 1, max))
  mae           = mean(abs(predicted - actual))
  mape          = mean(abs(predicted - actual)/actual)
  mse           = mean((actual - predicted)^2)
  rmse          = sqrt(mse)
  nrmse_mean    = rmse/mean(actual)
  cv_prcnt      = nrmse_mean*100 
  Var           = sum((actual - mean(actual))^2)
  RSS           = sum((actual - predicted)^2)
  var_r_squared = 1 - RSS / Var
  cv_prcnt      = nrmse_mean*100
     
     if(statistic=="EfronRSquared")
      {STAT=signif(var_r_squared, digits=digits)
      names(STAT)="EfronRSquared"}
     
     if(statistic=="MinMaxAccuracy")
      {STAT=signif(mma, digits=digits)
      names(STAT)="MinMaxAccuracy"}
     
     if(statistic=="MAE")
      {STAT=signif(mae, digits=digits)
      names(STAT)="MAE"}
     
     if(statistic=="MAPE")
      {STAT=signif(mape, digits=digits)
      names(STAT)="MAPE"}
     
     if(statistic=="MSE")
      {STAT=signif(mse, digits=digits)
      names(STAT)="MSE"}
     
     if(statistic=="RMSE")
      {STAT=signif(rmse, digits=digits)
      names(STAT)="RMSE"}
  
     if(statistic=="NRMSE.Mean")
      {STAT=signif(nrmse_mean, digits=digits)
      names(STAT)="NRMSE.Mean"}   
  
     if(statistic=="CV")
      {STAT=signif(cv_prcnt, digits=digits)
      names(STAT)="CV"}
     
     if(plotit){plot(data$actual, data$predicted,
                      xlab="Actual", ylab="Predicted", ...)}
   
 return(STAT)       
 }