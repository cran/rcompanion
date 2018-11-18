#' @title Pseudo r-squared measures for various models
#'
#' @description Produces McFadden, Cox and Snell, and Nagelkerke pseudo 
#'              R-squared measures, along with p-values, for models.
#' 
#' @param fit The fitted model object for which to determine pseudo r-squared.
#' @param null The null model object against which to compare the fitted model 
#'             object. The null model must be nested in the fitted model to be 
#'             valid. Specifying the null
#'             is optional for some model object types
#'             and is required for others.
#' @param restrictNobs If \code{TRUE}, limits the observations for the null
#'                     model to those used in the fitted model.
#'                     Works with only some model object types.
#' @details  Pseudo R-squared values are not directly comparable to the  
#'           R-squared for OLS models.  Nor can they be interpreted as the  
#'           proportion of the variability in the dependent variable that is  
#'           explained by model. Instead pseudo R-squared measures are relative
#'           measures among similar models indicating how well the model
#'           explains the data.
#'           
#'           Cox and Snell is also referred to as ML. Nagelkerke is also  
#'           referred to as Cragg and Uhler.
#'           
#'           Model objects accepted are lm, glm, gls, lme, lmer, lmerTest, nls,
#'           clm, clmm, vglm, glmer, negbin, zeroinfl, betareg, and rq.            
#'                                       
#'           Model objects that require the null model to 
#'           be defined are nls, lmer, glmer, and clmm. 
#'           Other objects use the \code{update} function to
#'           define the null model.
#'           
#'           Likelihoods are found using ML (\code{REML = FALSE}).
#'           
#'           The fitted model and the null model
#'           should be properly nested.
#'           That is, the terms of one need to be a subset of the the other,
#'           and they should have the same set of observations.
#'           One issue arises when there are \code{NA}
#'           values in one variable but not another, and observations with 
#'           \code{NA} are removed in the model fitting.  The result may be
#'           fitted and null models with
#'           different sets of observations.
#'           Setting \code{restrictNobs} to \code{TRUE} 
#'           ensures that only observations in
#'           the fit model are used in the null model.
#'           This appears to work for \code{lm} and some \code{glm} models,
#'           but causes the function to fail for other model
#'           object types.
#'           
#'           Some pseudo R-squared measures may not be appropriate
#'           or useful for some model types.
#'           
#'           Calculations are based on log likelihood values for models.
#'           Results may be different than those based on deviance.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/G_10.html}
#' @seealso \code{\link{accuracy}}
#' @concept pseudo r-squared cox snell nagelkerke likelihood
#' @return A list of six objects describing the models used, the pseudo 
#'         r-squared values, the likelihood ratio test for the model,
#'         the number of obervaton for the models,
#'         messages, and any warnings.
#'
#' @section Acknowledgements:
#'          My thanks to
#'          Jan-Herman Kuiper of Keele University for suggesting
#'          the \code{restrictNobs} fix.
#'          
#' @examples
#' ### Logistic regression example
#' data(AndersonBias)
#' model = glm(Result ~ County + Sex + County:Sex,
#'            weight = Count,
#'            data = AndersonBias,
#'            family = binomial(link="logit"))
#' nagelkerke(model)
#' 
#' ### Quadratic plateau example 
#' ### With nls, the  null needs to be defined
#' data(BrendonSmall)
#' quadplat = function(x, a, b, clx) {
#'           ifelse(x  < clx, a + b * x   + (-0.5*b/clx) * x   * x,
#'                            a + b * clx + (-0.5*b/clx) * clx * clx)}
#' model = nls(Sodium ~ quadplat(Calories, a, b, clx),
#'             data = BrendonSmall,
#'             start = list(a   = 519,
#'                          b   = 0.359,
#'                          clx = 2304))
#' nullfunct = function(x, m){m}
#' null.model = nls(Sodium ~ nullfunct(Calories, m),
#'              data = BrendonSmall,
#'              start = list(m   = 1346))
#' nagelkerke(model, null=null.model)
#' 
#' @importFrom stats df.residual logLik nobs pchisq update nls
#' 
#' @export

nagelkerke = 
function(fit, null=NULL, restrictNobs=FALSE)
{
   TOGGLE =   (class(fit)[1]=="lm"
             | class(fit)[1]=="gls"
             | class(fit)[1]=="lme"
             | class(fit)[1]=="glm"
             | class(fit)[1]=="negbin"
             | class(fit)[1]=="zeroinfl"
             | class(fit)[1]=="clm"
             | class(fit)[1]=="vglm"
             | class(fit)[1]=="betareg"
             | class(fit)[1]=="rq")
   BOGGLE =   (class(fit)[1]=="nls"
             | class(fit)[1]=="lmerMod"
             | class(fit)[1]=="glmerMod"
             | class(fit)[1]=="merModLmerTest"
             | class(fit)[1]=="lmerModLmerTest"
             | class(fit)[1]=="clmm")
   SMOGGLE =   (class(fit)[1]=="lmerMod"
              | class(fit)[1]=="glmerMod"
              | class(fit)[1]=="merModLmerTest"
              | class(fit)[1]=="lmerModLmerTest"
              | class(fit)[1]=="vglm")
   ZOGGLE  = (class(fit)[1]=="zeroinfl")
   ZOGGLE2 = (class(fit)[1]=="rq")
   NOGGLE = is.null(null)
   ERROR  = "Note: For models fit with REML, these statistics are based on refitting with ML"
   ERROR2 = "None"
   
  if(!restrictNobs & NOGGLE  & TOGGLE){null = update(fit, ~ 1)}
  if(restrictNobs  & NOGGLE  & TOGGLE){null = update(fit, ~ 1, data=fit$model)}
   
  if(restrictNobs  & !NOGGLE){null = update(null, data=fit$model)}
    
  if(NOGGLE & BOGGLE)
     {ERROR = "You need to supply a null model for nls, lmer, glmer, or clmm"}
  if((!TOGGLE) & (!BOGGLE))
   {ERROR = "This function will work with lm, gls, lme, lmer, glmer, glm, negbin, zeroinfl, nls, clm, clmm, and vglm"}
  
   SMOGGLE2 = (class(null)[1]=="lmerMod"
              | class(null)[1]=="glmerMod"
              | class(null)[1]=="merModLmerTest"
              | class(null)[1]=="lmerModLmerTest"
              | class(null)[1]=="vglm")   
   
  Y = matrix(rep(NA,2),
            ncol=1)
  colnames(Y) = ""
  rownames(Y) = c("Model:", "Null:")
  
  Z = matrix(rep(NA, 3),
             ncol=1)
  colnames(Z) = c("Pseudo.R.squared")
  rownames(Z) = c("McFadden", "Cox and Snell (ML)", 
                  "Nagelkerke (Cragg and Uhler)")
  
  X = matrix(rep(NA,4),
             ncol=4)
  colnames(X) = c("Df.diff","LogLik.diff","Chisq","p.value")
  rownames(X) = ""
  
  U = matrix(rep(NA,2),
            ncol=1)
  colnames(U) = ""
  rownames(U) = c("Model:", "Null:")
  
  if(TOGGLE | BOGGLE){
  if (!SMOGGLE){Y[1]= toString(fit$call)}
  if (SMOGGLE){Y[1]= toString(fit@call)}
  }
 
  if(TOGGLE | (BOGGLE & !NOGGLE)){
 
  if (!SMOGGLE2){Y[2]= toString(null$call)}
  if (SMOGGLE2){Y[2]= toString(null@call)}
 
  if(!ZOGGLE & !ZOGGLE2){N = nobs(fit)
                         U[1,1]= nobs(fit); U[2,1]= nobs(null)}
  if(!ZOGGLE &  ZOGGLE2){N = length(fit$y)
                         U[1,1]= length(fit$y); U[2,1]= length(null$y)}
  if(ZOGGLE){N = fit$n
             U[1,1]= fit$n; U[2,1]= null$n}

  if(U[1,1] != U[2,1]){
    ERROR2 = "WARNING: Fitted and null models have different numbers of observations"}
  
  m = suppressWarnings(logLik(fit, REML=FALSE))[1]
  n = suppressWarnings(logLik(null, REML=FALSE))[1]
  mf = 1 - m/n
  Z[1,] = signif(mf, digits=6)
  cs = 1 - exp(-2/N * (m - n))
  Z[2,] = signif(cs, digits=6)
  nk = cs/(1 - exp(2/N * n))
  Z[3,] = signif(nk, digits=6)
  
  o = n - m
  dfm = attr(logLik(fit),"df")
  dfn = attr(logLik(null),"df")
  if(class(fit)[1]=="vglm"){dfm=df.residual(fit)}
  if(class(fit)[1]=="vglm"){dfn=df.residual(null)}
  dff = dfn - dfm
  CHI = 2 * (m - n)
  P = pchisq(CHI, abs(dff), lower.tail = FALSE)
  
  X [1,1] = dff
  X [1,2] = signif(o, digits=5)             
  X [1,3] = signif(CHI, digits=5)
  X [1,4] = signif(P, digits=5)     
  }
  
  W=ERROR
  
  WW=ERROR2
  
  V = list(Y, Z, X, U, W, WW) 
  names(V) = c("Models", "Pseudo.R.squared.for.model.vs.null", 
               "Likelihood.ratio.test", "Number.of.observations",
               "Messages", "Warnings")
  return(V)            
}