#' @title Compare fit statistics for lm models
#'
#' @description Produces a table of fit statistics for multiple lm models.
#' 
#' @param fits A series of model object names, separated by commas.
#' @param ... Other arguments passed to \code{list}.
#' 
#' @details  Produces a table of fit statistics for multiple lm models: 
#'           AIC, AICc, BIC, p-value, R-squared, and adjusted R-squared.
#'           
#'           Smaller values for AIC, AICc, and BIC indicate a better balance
#'           of goodness-of-fit of the model and the complexity of the
#'           model. The goal is to find a model that adequately explains the 
#'           data without having too many terms.
#'           
#'           BIC tends to choose models with fewer parameters relative to AIC.
#'           
#'           In the table, \code{Shapiro.W} and \code{Shapiro.p} are the
#'           W statistic and p-value for the Shapiro-Wilks test on the residuals
#'           of the model.
#'           
#'           For comparisons with AIC, etc., to be valid, both models must
#'           have the same data, without transformations, use the same 
#'           dependent variable, and be fit with the same method.
#'           They do not need to be nested.
#'           
#'           The function will fail if a model formula is
#'           longer than 500 characters.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/I_10.html}
#'             \url{http://rcompanion.org/rcompanion/e_05.html}
#' @seealso \code{\link{compareGLM}}, \code{\link{pairwiseModelAnova}}, , \code{\link{accuracy}}
#' @concept AIC BIC AICc p-value R-squared
#' @return A list of two objects: The series of model calls, and a data 
#'         frame of statistics for each model.
#'         
#' @examples
#' ### Compare among polynomial models
#' data(BrendonSmall)
#' BrendonSmall$Calories = as.numeric(BrendonSmall$Calories)
#'
#' BrendonSmall$Calories2 = BrendonSmall$Calories * BrendonSmall$Calories
#' BrendonSmall$Calories3 = BrendonSmall$Calories * BrendonSmall$Calories * 
#'                          BrendonSmall$Calories
#' BrendonSmall$Calories4 = BrendonSmall$Calories * BrendonSmall$Calories * 
#'                          BrendonSmall$Calories * BrendonSmall$Calories
#' model.1 = lm(Sodium ~ Calories, data = BrendonSmall)
#' model.2 = lm(Sodium ~ Calories + Calories2, data = BrendonSmall)
#' model.3 = lm(Sodium ~ Calories + Calories2 + Calories3, data = BrendonSmall)
#' model.4 = lm(Sodium ~ Calories + Calories2 + Calories3 + Calories4,
#'              data = BrendonSmall)
#' compareLM(model.1, model.2, model.3, model.4)
#' 
#' @importFrom stats logLik nobs pf shapiro.test
#' 
#' @export

compareLM = 
function (fits, ...) 
{
 fits = list(fits, ...)
 n = length(fits)
 Y = matrix(rep("A",n),
            ncol=1)
 colnames(Y) = "Formula"
 rownames(Y) = seq(1,n)

  for(i in 1:n)
    {
     Y[i,]= deparse(formula(fits[[i]]), width.cutoff = 500L)
     }
   
 Z = data.frame(Rank=rep(NA,n),
                Df.res=rep(NA,n),
                AIC=rep(NA,n),
                AICc = rep(NA,n),
                BIC=rep(NA,n),
                R.squared=rep(NA,n),
                Adj.R.sq=rep(NA,n),
                p.value=rep(NA,n),
                Shapiro.W=rep(NA,n),
                Shapiro.p=rep(NA,n),
                stringsAsFactors=FALSE)   
   for(i in 1:n)
    {
     k = length(fits[[i]]$coefficients)+1
     L = as.numeric(logLik(fits[[i]]))
     N = nobs(fits[[i]])
     aic  = -2*L+2*k
     aicc = aic+2*(k*(k+1))/(N-k-1)
     bic  = -2*L+log(N)*k
    
     Z[i,]=c(signif(fits[[i]]$rank, digits=4),
             signif(fits[[i]]$df.residual, digits=4),
             signif(aic, digits=4),
             signif(aicc, digits=4),
             signif(bic, digits=4),
             signif(summary(fits[[i]])$r.squared, digits=4),
             signif(summary(fits[[i]])$adj.r.squared, digits=4),
             signif(pf(summary(fits[[i]])$fstatistic[1], 
                summary(fits[[i]])$fstatistic[2], 
                summary(fits[[i]])$fstatistic[3],
                lower.tail = FALSE), digits=4),
             signif(as.numeric(shapiro.test(fits[[i]]$residuals)$statistic),
                    digits =4),
             signif(shapiro.test(fits[[i]]$residuals)$p.value,
                    digits=4)
             )
     }     
   
 W = list(Y, Z)
 names(W) = c("Models",
              "Fit.criteria")
 return(W)            
 } 