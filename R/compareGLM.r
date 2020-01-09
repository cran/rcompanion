#' @title Compare fit statistics for glm models
#'
#' @description Produces a table of fit statistics for multiple glm models.
#' 
#' @param fits A series of model object names, separated by commas.
#' @param ... Other arguments passed to \code{list}.
#' 
#' @details  Produces a table of fit statistics for multiple glm models: 
#'           AIC, AICc, BIC, p-value, pseudo R-squared
#'           (McFadden, Cox and Snell, Nagelkerke).
#'           
#'           Smaller values for AIC, AICc, and BIC indicate a better balance
#'           of goodness-of-fit of the model and the complexity of the
#'           model. The goal is to find a model that adequately explains the 
#'           data without having too many terms.
#'           
#'           BIC tends to choose models with fewer parameters relative to AIC.

#'           For comparisons with AIC, etc., to be valid, both models must
#'           have the same data, without transformations, use the same 
#'           dependent variable, and be fit with the same method.
#'           They do not need to be nested.
#'           
#'           The function will fail if a model formula is
#'           longer than 500 characters.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/rcompanion/e_07.html}
#' @seealso \code{\link{compareLM}}, \code{\link{pairwiseModelAnova}}, \code{\link{accuracy}}
#' @concept AIC BIC AICc p-value pseudo R-squared glm
#' @return A list of two objects: The series of model calls, and a data 
#'         frame of statistics for each model.
#'         
#' @examples
#' ### Compare among logistic regresion models
#' data(AndersonBias)
#' model.0 = glm(Result ~ 1, weight = Count, data = AndersonBias,
#'              family = binomial(link="logit"))
#' model.1 = glm(Result ~ County, weight = Count, data = AndersonBias,
#'              family = binomial(link="logit"))
#' model.2 = glm(Result ~ County + Sex, weight = Count, data = AndersonBias,
#'              family = binomial(link="logit"))
#' model.3 = glm(Result ~ County + Sex + County:Sex, weight = Count, 
#'              data = AndersonBias, family = binomial(link="logit"))
#' compareGLM(model.0, model.1, model.2, model.3)              
#' 
#' @importFrom stats logLik nobs dchisq update
#' 
#' @export
 
compareGLM = 
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
                McFadden=rep(NA,n),
                Cox.and.Snell=rep(NA,n),
                Nagelkerke=rep(NA,n),
                p.value=rep(NA,n),
                stringsAsFactors=FALSE)   
   for(i in 1:n)
    {
     k    = length(fits[[i]]$coefficients)+1
     L    = as.numeric(logLik(fits[[i]]))
     N    = nobs(fits[[i]])
     m    = logLik(fits[[i]])[1]
     null = update(fits[[i]], ~1)
     n    = logLik(null)[1]
     aic  = -2*L+2*k
     aicc = aic+2*(k*(k+1))/(N-k-1)
     bic  = -2*L+log(N)*k
     mf   = 1 - m/n
     cs   = 1 - exp(-2/N * (m - n))
     nk   = cs/(1 - exp(2/N * n))
     p.val= dchisq(fits[[i]]$null.deviance-fits[[i]]$deviance, 
            fits[[i]]$df.null-fits[[i]]$df.residual)

     
     Z[i,]=c(signif(fits[[i]]$rank, digits=4),
             signif(fits[[i]]$df.residual, digits=4),
             signif(aic, digits=4),
             signif(aicc, digits=4),
             signif(bic, digits=4),
             signif(mf, digits=4),
             signif(cs, digits=4),
             signif(nk, digits=4),
             signif(p.val, digits=4)
             )
     }     
   
 W = list(Y, Z)
 names(W) = c("Models",
              "Fit.criteria")
 return(W)          
 }