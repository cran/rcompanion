#' @title Pairwise permutation tests for percentiles
#'
#' @description Conducts pairwise permutation tests across groups 
#'              for percentiles, medians, and proportion below a
#'              threshold value.
#' 
#' @param formula A formula indicating the response variable and
#'                the independent variable. e.g. y ~ group.
#' @param data   The data frame to use. 
#' @param x If no formula is given, the response variable for one group.
#' @param y The response variable for the other group.
#' @param test The statistic to compare between groups.  Can be
#'             \code{"median"}, \code{"percentile"}, \code{"iqr"},
#'             \code{"proportion"},
#'             \code{"mean"}, or \code{"variance"}.  
#' @param tau If \code{"percentile"} is chosen as the \code{test},
#'            \code{tau} indicates the percentile to test.  Expressed
#'            as a quantile.  That is, 0.5 indicates a test for medians.
#'            0.75 indicates a test for 75th percentiles.
#' @param type The \code{type} value passed to the \code{quantile} function.
#' @param threshold If \code{"proportion"} is chosen as the \code{test},
#'            \code{threshold} indicates the value of the dependent variable
#'            to use as the threshold.  For example, to test if there is a 
#'            different in the proportion of observations below $10,000,
#'            \code{threshold = 10000} would be used.
#' @param comparison If \code{"proportion"} is chosen as the \code{test},
#'            \code{comparison} indicates the inequality to use.  Options are
#'            \code{"<"}, \code{"<="}, \code{">"}, \code{">="}, or , \code{"=="}
#' @param r The number of replicates in the permutation test.
#' @param method The p-value adjustment method to use for multiple tests.
#'               See \code{\link{p.adjust}}.
#' @param digits The number of significant digits in the output. 
#' @param progress If \code{TRUE}, prints a dot for every 1 percent of the 
#'                 progress while conducting the test.
#'             
#' @details The function conducts pairwise tests using the 
#'          \code{percentileTest} function. The user can consult the
#'          documentation for that function for additional details.
#'          
#'          The input should include either \code{formula} and \code{data};
#'          or \code{x}, and \code{y}.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_15.html}
#' @seealso \code{\link{percentileTest}}, \code{\link{groupwisePercentile}}
#' @concept median percentile quantile permutation
#' @return A dataframe of the groups being compared, the p-values,
#'         and the adjusted p-values.
#'         
#' @note    The parsing of the formula is simplistic. 
#'          The first variable on the
#'          left side is used as the measurement variable.  
#'          The first variable on the
#'          right side is used for the grouping variable.   
#'          
#' @examples
#' data(BrendonSmall)
#' PT = pairwisePercentileTest(Sodium ~ Instructor, 
#'                             data = BrendonSmall, 
#'                             test = "percentile", 
#'                             tau  = 0.75)
#' PT
#' cldList(p.adjust ~ Comparison,
#'         data       = PT,
#'         threshold  = 0.05)
#'         
#' data(BrendonSmall)
#' PT = pairwisePercentileTest(Sodium ~ Instructor, 
#'                             data       = BrendonSmall, 
#'                             test       = "proportion", 
#'                             threshold  = 1300)
#' PT
#' cldList(p.adjust ~ Comparison,
#'         data       = PT,
#'         threshold  = 0.05)                         
#' 
#' @importFrom stats p.adjust
#' 
#' @export

pairwisePercentileTest = 
 function(formula=NULL, data=NULL,
           x=NULL, y=NULL,
           test="median", tau=0.5, type=7,
           threshold = NA, comparison = "<",
           r=1000, digits=4, progress="TRUE",
           method = "fdr")
  {
  if(!is.null(formula)){
    x  = eval(parse(text=paste0("data","$",all.vars(formula[[2]])[1])))
    g  = eval(parse(text=paste0("data","$",all.vars(formula[[3]])[1])))
  }
  if(is.factor(g)){g=droplevels(g)}
  if(!is.factor(g)){g=factor(g)}
  n = length(levels(g))
  N = n*(n-1)/2
  d = data.frame(x = x, g = g)
  Z = data.frame(Comparison=rep("A", N),
                 p.value=rep(NA, N),
                 p.adjust=rep(NA, N),
                 stringsAsFactors=FALSE)
  k=0               
  for(i in 1:(n-1)){
     for(j in (i+1):n){
       k=k+1
    if(progress){cat("Comparison ", k, "\n")}
     Namea = as.character(levels(g)[i])
     Nameb = as.character(levels(g)[j])
     Datax = subset(d, g==levels(g)[i])
     Datay = subset(d, g==levels(g)[j])
     Dataz = rbind(Datax, Datay)
     Dataz$g2 = droplevels(Dataz$g)
     
     z = percentileTest(x ~ g2, data=Dataz,
                        test=test, tau=tau, 
                        type=type,
                        threshold = threshold, 
                        comparison = comparison,
                        r=r, digits=digits, 
                        progress=progress)
    
     P = signif(z[["Result"]]$p.value, digits=digits)
     P.adjust = NA                       
     Z[k,] =c( paste0(Namea, " - ", Nameb, " = 0"), 
             P, P.adjust)
     }
    } 
  Z$p.adjust = signif(p.adjust(Z$p.value, method = method), digits=digits) 
  return(Z)
  }