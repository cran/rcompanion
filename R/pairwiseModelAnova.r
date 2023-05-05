 #' @title Compare model objects with F test and likelihood ratio test
#'
#' @description Compares a series of models with pairwise F tests and
#'              likelihood ratio tests.
#' 
#' @param fits A series of model object names, separated by commas.
#' @param ... Other arguments passed to \code{list}.
#' 
#' @details  For comparisons to be valid, both models must
#'           have the same data, without transformations, use the same 
#'           dependent variable, and be fit with the same method.
#'           
#'           To be valid, models need to be nested.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @seealso \code{\link{compareGLM}}, \code{\link{compareLM}}
#' 
#' @concept likelihood ratio test
#' @concept F test
#' 
#' @return A list of: The calls of the models compared; a data frame of 
#'         comparisons and F tests; and a data frame of
#'         comparisons and likelihood ratio tests.
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
#' pairwiseModelAnova(model.1, model.2, model.3, model.4)
#' 
#' @importFrom stats logLik 
#' @importFrom lmtest lrtest
#' 
#' @export

pairwiseModelAnova = 
  function (fits, ...)
  {
  fits = list(fits, ...)
  n = length(fits)
  N = n*(n-1)/2
  
  Y = data.frame(Model=rep("A",n),
                 Call=rep("A",n),
                 stringsAsFactors=FALSE)
for(i in 1:n)
    {
     Y[i,]= deparse(formula(fits[[i]]))
     }
  
  Z = data.frame(Comparison = rep("A", N),
                 Df.diff    = rep(NA, N),
                 RSS.diff   = rep(NA, N),
                 F          = rep(NA, N),
                 p.value    = rep(NA, N),
                 stringsAsFactors=FALSE)
  k=0              
  for(i in 1:(n-1)){
     for(j in (i+1):n){
       k=k+1
       Namea = as.character(i)
       Nameb = as.character(j) 
       av = anova(fits[[i]], fits[[j]])
          
       Z[k,] =c(paste0(Namea, " - ", Nameb, " = 0"), 
                as.numeric(av$Df[2]),
                signif(as.numeric(av$"Sum of Sq"[2], digits=4)),
                signif(as.numeric(av$F[2], digits=4)),
                signif(as.numeric(av$"Pr(>F)"[2], digits=4))
                )
       }
     }
             
  X = data.frame(Comparison = rep("A", N),
                 Df.diff    = rep(NA, N),
                 LogLik.diff= rep(NA, N),
                 Chi.sq     = rep(NA, N),
                 p.value    = rep(NA, N),
                 stringsAsFactors=FALSE)
  k=0              
  for(i in 1:(n-1)){
     for(j in (i+1):n){
       k=k+1
       Namea = as.character(i)
       Nameb = as.character(j) 
       lr = lrtest(fits[[i]], fits[[j]])
          
       X[k,] =c(paste0(Namea, " - ", Nameb, " = 0"), 
                as.numeric(lr$Df[2]),
                signif(as.numeric(lr$LogLik[2])-as.numeric(lr$LogLik[1]), 
                       digits=4),
                signif(as.numeric(lr$Chisq[2]), digits=4),
                signif(as.numeric(lr$"Pr(>Chisq)"[2]), digits=4)
                )
       }
     }
              
 W = list(Y, Z, X)
 names(W) = c("Models",
              "Anova.test.table",
              "Likelihood.ratio.test.table")
 return(W)            
 }
     
  