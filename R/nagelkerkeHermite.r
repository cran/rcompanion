#' @title Pseudo r-squared measures for hermite models
#'
#' @description Produces McFadden, Cox and Snell, and Nagelkerke pseudo 
#'              R-squared measures, along with p-value for the model, 
#'              for hermite regression objects.
#' 
#' @param fit The fitted model object for which to determine pseudo r-squared.
#' @param null The null model object against which to compare the fitted model 
#'             object. The null model must be nested in the fitted model to be 
#'             valid.
#' @details  Hermite regression is performed with the \code{hermite} package.
#'           
#'           For pseudo r-squared measures, Cox and Snell is also referred to 
#'           as ML. Nagelkerke is also referred to as Cragg and Uhler.
#'           
#'           The fit model and the null model
#'           should be properly nested.
#'           That is, the terms of one need to be a subset of the the other,
#'           and they should have the same set of observations.           
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/J_01.html}
#' @seealso \code{\link{nagelkerke}}
#' @concept pseudo r-squared cox snell nagelkerke likelihood hermite
#' @return A list of six objects describing the models used, the pseudo 
#'         r-squared values, the likelihood ratio test for the model, the AIC
#'         for the fitted and null models,
#'         the number of observations for the models,
#'         and any warnings.
#'         
#' @examples
#' data(Monarchs)
#' library(hermite)
#' model = glm.hermite(Monarchs ~ Garden,
#'                     data = Monarchs,
#'                     link = "log",
#'                     m=3)
#' null = glm.hermite(Monarchs ~ 1,
#'                    data = Monarchs,
#'                    link = "log",
#'                    m=3)
#' nagelkerkeHermite(model, null)
#' 
#' @importFrom hermite glm.hermite
#' @importFrom stats df.residual logLik nobs pchisq update
#' 
#' @export

nagelkerkeHermite = 
function(fit, null)
{  
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
  
  Y = matrix(rep("MISSING",2),
            ncol=1)
  colnames(Y) = ""
  rownames(Y) = c("Model:", "Null:")
  
  U = matrix(rep(NA,2),
            ncol=1)
  colnames(U) = ""
  rownames(U) = c("Model:", "Null:")
  
  WW = "None"
  
  Y[1] = toString(summary(fit)$call)
  Y[2] = toString(summary(null)$call)
  U[1] = length(fit$fitted.values)
  U[2] = length(null$fitted.values)
  
  if (U[1] != U[2]){
    WW = "WARNING: Fitted and null models have different numbers of observations"}
 
  N = length(fit$fitted.values)
  m = fit$loglik
  n = null$loglik
  mf = 1 - m/n
  Z[1,] = signif(mf, digits=6)
  cs = 1 - exp(-2/N * (m - n))
  Z[2,] = signif(cs, digits=6)
  nk = cs/(1 - exp(2/N * n))
  Z[3,] = signif(nk, digits=6)
 
  o = n - m
  dff =  length(null$coefs)- length(fit$coefs)
  CHI = 2 * (m - n)
  P = pchisq(CHI, abs(dff), lower.tail = FALSE)
  
  X [1,1] = dff
  X [1,2] = signif(o, digits=5)             
  X [1,3] = signif(CHI, digits=5)
  X [1,4] = signif(P, digits=5)     
  
  W = matrix(rep(NA,2),
             ncol=1)
  colnames(W) = "AIC"
  rownames(W) = c("Model:", "Null:")
  
  W[1] = signif(summary(fit)$aic, digits=5)
  W[2] = signif(summary(null)$aic, digits=5)
  
  V = list(Y, Z, X, W, U, WW) 
  names(V) = c("Models", "Pseudo.R.squared.for.model.vs.null", "Likelihood.ratio.test", "AIC", "Number.of.observations",
               "Warnings")
  return(V)            
}