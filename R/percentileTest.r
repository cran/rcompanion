#' @title Test of percentiles by permutation test
#'
#' @description Conducts a permutation test to compare two groups for medians,
#'              percentiles, or proportion below a threshold value.
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
#' @param digits The number of significant digits in the output. 
#' @param progress If \code{TRUE}, prints a dot for every 1 percent of 
#'                 progress while conducting the test.
#' 
#' @details The function will test for a difference in medians, percentiles,
#'          interquartile ranges, proportion of observations above or below
#'          some threshold value, 
#'          means, or variances between two groups
#'          by permutation test.
#'          
#'          The permutation test simply permutes the observed values over the
#'          two groups and counts how often the calculated statistic is
#'          at least as extreme as the original observed statistic.
#'          
#'          The input should include either \code{formula} and \code{data};
#'          or \code{x} and \code{y}.
#'          
#'          The function removes cases with NA in any of the variables.
#'          
#'          If the independent variable has more than two groups,
#'          only the first two levels of the factor variable will be used.
#'          
#'          The p-value returned is a two-sided test.
#'           
#' @note    The parsing of the formula is simplistic. 
#'          The first variable on the
#'          left side is used as the measurement variable.  
#'          The first variable on the
#'          right side is used for the independent variable.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{http://rcompanion.org/handbook/F_15.html}
#'             
#' @concept percentile
#' @concept median
#' @concept quantile
#' @concept permutation test
#' 
#' @return A list of three data frames with
#'         the data used, a summary for each group,
#'         and the p-value from the test.
#'  
#' @examples
#' data(BrendonSmall)
#' percentileTest(Sodium ~ Instructor, 
#'                data=BrendonSmall, 
#'                test="median")
#'
#' percentileTest(Sodium ~ Instructor, 
#'                data=BrendonSmall, 
#'                test="percentile", 
#'                tau = 0.75)
#'
#' percentileTest(Sodium ~ Instructor, 
#'                data=BrendonSmall, 
#'                test="proportion", 
#'                threshold = 1300)
#'                
#' @importFrom utils tail
#' @importFrom stats median var quantile IQR complete.cases
#' 
#' @export

percentileTest = 
  function(formula=NULL, data=NULL,
           x=NULL, y=NULL,
           test="median", tau=0.5, type=7,
          threshold = NA, comparison = "<",
           r=1000, digits=4, progress="TRUE")
  {
  if(is.null(formula)){
    xy   = c(x,y)
    xname = paste(as.character(substitute(x)), collapse=" ")
    yname = paste(as.character(substitute(y)), collapse=" ")
    g     = factor(c(rep(xname, length(x)), rep(yname, length(y))))
    }
if(!is.null(formula)){
xy   = eval(parse(text=paste0("data","$",all.vars(formula[[2]])[1])))
g    = eval(parse(text=paste0("data","$",all.vars(formula[[3]])[1])))
if(!is.factor(g)){g=factor(g)}
}
 Complete = complete.cases(xy, g)
 xy = xy[Complete]
 g   = g[Complete]
 g   = droplevels(g)
 x   = xy[g == levels(g)[1]]
 y   = xy[g == levels(g)[2]]
 xy  = c(x, y)
 
 PropA = NA
 PropB = NA
 
 if(test=="mean"){      Diff = abs(mean(x) - mean(y))}
 if(test=="median"){    Diff = abs(median(x) - median(y))}
 if(test=="variance"){  Diff = abs(var(x) - var(y))}
 if(test=="percentile"){Diff = abs(quantile(x, probs=tau, type=type) -
                                   quantile(y, probs=tau, type=type))}
 if(test=="iqr"){       Diff = abs(IQR(x, type=type) -
                                   IQR(y, type=type))}
 if(test=="proportion"){
      if(is.na(threshold)){threshold = median(c(x,y))}
      SumA = sum(eval(parse(text=paste0("x", comparison, threshold))))
      SumB = sum(eval(parse(text=paste0("y", comparison, threshold))))
      PropA = SumA / length(x)
      PropB = SumB / length(y)
      Diff = abs(PropA - PropB)
      }
  Count = 0
  if (progress){tick=r/100}
  for(i in 1:r){
  S = sample(xy, replace=FALSE)
  Sx   = S[g == levels(g)[1]]
  Sy   = S[g == levels(g)[2]]

  if(test=="mean"){      Sdiff = abs(mean(Sx) - mean(Sy))}
  if(test=="median"){    Sdiff = abs(median(Sx) - median(Sy))}
  if(test=="variance"){  Sdiff = abs(var(Sx) - var(Sy))}
  if(test=="percentile"){Sdiff = abs(quantile(Sx, probs=tau, type=type) -
                                     quantile(Sy, probs=tau, type=type))}
  if(test=="iqr"){       Sdiff = abs(IQR(Sx, type=type) -
                                     IQR(Sy, type=type))}
  if(test=="proportion"){
      SA = sum(eval(parse(text=paste0("Sx", comparison, threshold))))
      SB = sum(eval(parse(text=paste0("Sy", comparison, threshold))))
      PA = SA / length(Sx)
      PB = SB / length(Sy)
      Sdiff = abs(PA -PB)
      }
  if(Sdiff >= Diff){Count=Count+1}
  if(progress){if(i%%tick==0){cat(".")}}
  }
  if(progress){cat("\n"); cat("\n")}
  P.value = Count/r
  Z = data.frame(
    Statistic  = rep("NA", 2),
    n          = rep(as.numeric(NA), 2),
    mean       = rep(as.numeric(NA), 2),
    sd         = rep(as.numeric(NA), 2),
    min        = rep(as.numeric(NA), 2),
    p25        = rep(as.numeric(NA), 2),
    median     = rep(as.numeric(NA), 2),
    p75        = rep(as.numeric(NA), 2),
    max        = rep(as.numeric(NA), 2),
    iqr        = rep(as.numeric(NA), 2),
    comparison = rep("NA", 2),
    threshold  = rep(as.numeric(NA), 2),
    proportion = rep(as.numeric(NA), 2),
    stringsAsFactors=FALSE)
  
  Z[1,1] = levels(g)[1]
  Z[1,2] = signif(length(x),digits)
  Z[1,3] = signif(mean(x),digits)
  Z[1,4] = signif(sd(x),digits)
  Z[1,5] = signif(min(x),digits)
  Z[1,6] = signif(quantile(x,0.25),digits)
  Z[1,7] = signif(median(x),digits)
  Z[1,8] = signif(quantile(x,0.75),digits)
  Z[1,9] = signif(max(x),digits)
  Z[1,10] = signif(IQR(x),digits)
  Z[1,11] = comparison
  Z[1,12] = signif(threshold, digits)
  Z[1,13] = signif(PropA, digits)
            
  Z[2,1] = levels(g)[2]
  Z[2,2] = signif(length(y),digits)
  Z[2,3] = signif(mean(y),digits)
  Z[2,4] = signif(sd(y),digits)
  Z[2,5] = signif(min(y),digits)
  Z[2,6] = signif(quantile(y,0.25),digits)
  Z[2,7] = signif(median(y),digits)
  Z[2,8] = signif(quantile(y,0.75),digits)
  Z[2,9] = signif(max(y),digits)
  Z[2,10] = signif(IQR(y),digits)
  Z[2,11] = comparison
  Z[2,12] = signif(threshold, digits)
  Z[2,13] = signif(PropB, digits)
  
  colnames(Z)[1]=""
  if(test!="proportion"){Z = Z[, !names(Z) %in% 
     c('comparison', 'threshold', 'proportion')]}
  
  U = data.frame(
      Statistic = rep("NA", 1),
      p.value   = rep(as.numeric(NA), 1),
      stringsAsFactors=FALSE)
  U[1,1] = "p-value"
  U[1,2] = signif(P.value, digits)
  colnames(U)[1]=""
  
  V = data.frame(
      Formula   = "NA",
      Data      = "NA",
      Test      = "NA",
      tau       = as.numeric(NA),
      stringsAsFactors=FALSE)

  if(!is.null(formula)){V[1,1] = deparse(formula)}
  if(!is.null(formula)){V[1,2] = as.character(substitute(data))}
  V[1,3] = test
  V[1,4] = tau
  if(is.null(formula)){V[1,1] = xname 
                       V[1,2] = yname
                       colnames(V)[1:2]=c("x", "y")
                       }
  
  if(test!="percentile"){V = V[, !names(V) %in% c('tau')]}
  
  W = list(Test=V,
           Summary=Z, 
           Result=U)
  
  return(W)
}