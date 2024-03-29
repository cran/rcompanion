#' @title Pairwise McNemar and related tests for Cochran Q test post-hoc
#'
#' @description Conducts pairwise McNemar, exact, and permutation
#'              tests as a post-hoc to Cochran Q test.
#'              
#' @param formula A formula indicating the measurement variable and
#'                the grouping variable. e.g. y ~ group | block.
#' @param data   The data frame to use. 
#' @param x The response variable.
#' @param g The grouping variable.
#' @param block The blocking variable.
#' @param test If \code{"exact"}, conducts an exact test of symmetry
#'             analogous to a McNemar test.
#'             If \code{"mcnemar"}, conducts a McNemar test of symmetry.
#'             If \code{"permutation"}, conducts a permutation test 
#'             analogous to a McNemar test.
#' @param method The method for adjusting multiple p-values.
#'               See \code{\link{p.adjust}}.
#' @param digits The number of significant digits in the output.
#' @param correct If \code{TRUE}, applies a continuity correction
#'                for the McNemar test.       
#' 
#' @details The component tables for the pairwise tests
#'          must be of size 2 x 2.
#'          
#'          The input should include either \code{formula} and \code{data};
#'          or \code{x}, \code{g}, and \code{block}.
#'      
#' @note    The parsing of the formula is simplistic. 
#'          The first variable on the
#'          left side is used as the measurement variable.  
#'          The first variable on the
#'          right side is used for the grouping variable.
#'          The second variable on the
#'          right side is used for the blocking variable.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{https://rcompanion.org/handbook/H_07.html}
#' 
#' @concept post-hoc
#' @concept Cochran's Q
#' @concept McNemar's test  
#' 
#' @return A list containing: a data frame of results of the global test;
#'         a data frame of results of the pairwise results;
#'         and a data frame mentioning the p-value adjustment method.

#'  
#' @seealso \code{\link{nominalSymmetryTest}}, 
#'          \code{\link{groupwiseCMH}},
#'          \code{\link{pairwiseNominalIndependence}}, 
#'          \code{\link{pairwiseNominalMatrix}}              
#'         
#' @examples
#' ### Cochran Q post-hoc example
#' data(HayleySmith)
#' library(DescTools)
#' CochranQTest(Response ~ Practice | Student,
#'              data = HayleySmith)
#' HayleySmith$Practice = factor(HayleySmith$Practice,
#'                           levels = c("MowHeight", "SoilTest",
#'                                      "Clippings", "Irrigation"))
#' PT = pairwiseMcnemar(Response ~ Practice | Student,
#'                      data    = HayleySmith,
#'                      test    = "exact",
#'                      method  = "fdr",
#'                      digits  = 3)
#' PT
#' PT = PT$Pairwise
#' cldList(comparison = PT$Comparison,
#'         p.value    = PT$p.adjust,
#'         threshold  = 0.05)
#'                                                              
#' @importFrom stats xtabs mcnemar.test binom.test
#' @importFrom DescTools CochranQTest
#' @importFrom coin symmetry_test pvalue
#' 
#' @export

pairwiseMcnemar = 
  function(formula=NULL, data=NULL,
           x=NULL, g=NULL, block=NULL, 
           test="exact", method="fdr", digits=3, correct=FALSE)
  {
  if(!is.null(formula)){
    x      = eval(parse(text=paste0("data","$",all.vars(formula[[2]])[1])))
    g      = eval(parse(text=paste0("data","$",all.vars(formula[[3]])[1])))
    block  = eval(parse(text=paste0("data","$",all.vars(formula[[3]])[2])))
  }
  if(!is.factor(g)){g=factor(g)}
  if(!is.factor(block)){block=factor(block)}  
  Name = as.character(levels(g))
  n = length(Name)
  N = n*(n-1)/2
  data = data.frame(x = x, g = g, block = block)
  
  if(test=="exact"){
  Z = data.frame(Comparison=rep("A", N),
                 Successes=rep(NA, N),
                 Trials=rep(NA, N),
                 p.value=rep(NA, N),
                 p.adjust=rep(NA, N),
                 stringsAsFactors=FALSE)
  }
  
  if(test=="mcnemar"){
  Z = data.frame(Comparison=rep("A", N),
                 chi.sq=rep(NA, N),
                 df=rep(NA, N),
                 p.value=rep(NA, N),
                 p.adjust=rep(NA, N),
                 stringsAsFactors=FALSE)
  }
  
  if(test=="permutation"){
  Z = data.frame(Comparison=rep("A", N),
                 Z=rep(NA, N),
                 p.value=rep(NA, N),
                 p.adjust=rep(NA, N),
                 stringsAsFactors=FALSE)
  }
  
  k=0
  for(i in 1:(n-1)){
     for(j in (i+1):n){
       k=k+1
     Namea = Name[i]
     Nameb = Name[j]
     Datax = subset(data, g==Name[i])
     Datay = subset(data, g==Name[j])
     Dataz = rbind(Datax, Datay)
     Dataz$g2 = factor(Dataz$g)
     
     Dataz$x = factor(Dataz$x)
     Dataz$x = as.numeric(Dataz$x)-1
     X = xtabs(x ~ g2 + block, data=Dataz)
     a = sum(X[1,]==1 & X[2,]==1)
     b = sum(X[1,]==1 & X[2,]==0)
     c = sum(X[1,]==0 & X[2,]==1)
     d = sum(X[1,]==0 & X[2,]==0)
     Y= matrix(c(a,b,c,d),ncol=2,byrow=TRUE)
     
     if(test=="exact"){
       Test=binom.test(b, (b+c), 0.5)
       P = signif(Test$p.value, digits=digits)
       Successes = Test$statistic
       Trials   = Test$parameter
       P.adjust = NA                       
       Z[k,] =c(paste0(Namea, " - ", Nameb, " = 0"), 
                Successes, Trials, P, P.adjust)
     }
     
     if(test=="mcnemar"){
       Test=mcnemar.test(Y, correct=correct)
       P = signif(Test$p.value, digits=digits)
       Chi = signif(Test$statistic, digits=digits)
       Df = Test$parameter
       P.adjust = NA                       
       Z[k,] =c(paste0(Namea, " - ", Nameb, " = 0"), 
                Chi, Df, P, P.adjust)
       }
     
if(test=="permutation"){
       Test=symmetry_test(as.table(Y))
       P = signif(pvalue(Test), digits=digits)
       z = as.numeric(signif(statistic(Test), digits=digits))
       P.adjust = NA          
       Z[k,] =c(paste0(Namea, " - ", Nameb, " = 0"), 
                z, P, P.adjust)
       }
     }
    } 
  Z$p.adjust = signif(p.adjust(Z$p.value, method = method), digits=digits) 
  
  U = data.frame(Test=test)
  V = data.frame(Method=method)
  W = list(U, V, Z)
  names(W) = c("Test.method",
               "Adustment.method",
               "Pairwise")
  return(W)
  }