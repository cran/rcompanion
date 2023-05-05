#' @title Dominance statistic for two-sample paired data
#'
#' @description Calculates a dominance
#'              effect size statistic
#'              for two-sample paired data
#'              with confidence intervals by bootstrap
#'
#' @param formula A formula indicating the response variable and
#'                the independent variable. e.g. y ~ group.
#' @param data   The data frame to use. 
#' @param x If no formula is given, the response variable for one group.
#' @param y The response variable for the other group.
#' @param ci If \code{TRUE}, returns confidence intervals by bootstrap.
#'           May be slow.
#' @param conf The level for the confidence interval.
#' @param type The type of confidence interval to use.
#'             Can be any of "\code{norm}", "\code{basic}", 
#'                           "\code{perc}", or "\code{bca}".
#'             Passed to \code{boot.ci}.
#' @param R The number of replications to use for bootstrap.
#' @param histogram If \code{TRUE}, produces a histogram of bootstrapped values.
#' @param digits The number of significant digits in the output.
#' @param na.rm If \code{TRUE}, removes \code{NA} values from
#'              the input vectors or data frame.
#' @param ... Additional arguments. 
#'             
#' @details The calculated \code{Dominance} statistic is simply 
#'          the proportion of observations in \code{x} greater the paired
#'          observations in \code{y},
#'          minus
#'          the proportion of observations in \code{x} less than the paired
#'          observations in \code{y}
#'          
#'          It will range from -1 to 1, with
#'          and 1 indicating that 
#'          the all the observations in \code{x} are greater than 
#'          the paired observations in \code{y},
#'          and -1 indicating that
#'          the all the observations in \code{y} are greater than 
#'          the paired observations in \code{x}.
#'          
#'          The input should include either \code{formula} and \code{data};
#'          or \code{x}, and \code{y}. If there are more than two groups,
#'          only the first two groups are used.
#'          
#'          This statistic is appropriate for truly ordinal data,
#'          and could be considered an effect size statistic for
#'          a two-sample paired sign test.
#'          
#'          Ordered category data need to re-coded as
#'          numeric, e.g. as with \code{as.numeric(Ordinal.variable)}.
#'
#'          When the statistic is close to 1 or close to -1,
#'          or with small sample size,
#'          the confidence intervals 
#'          determined by this
#'          method may not be reliable, or the procedure may fail.
#'           
#'          VDA is the analogous statistic, converted to a probability,
#'          ranging from 0 to 1, specifically,
#'          \code{VDA = Dominance / 2 + 0.5}
#'                      
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{https://rcompanion.org/handbook/F_07.html}
#' 
#' @seealso \code{\link{oneSampleDominance}},
#'          \code{\link{vda}},
#'          \code{\link{cliffDelta}}
#' 
#' @concept effect size
#' @concept dominance
#' @concept sign test
#' 
#' @return A small data frame consisting of descriptive statistics,
#'         the dominance statistic,
#'         and potentially the lower and upper confidence limits.
#'                            
#' @examples
#' data(Pooh)
#' Time.1 = Pooh$Likert[Pooh$Time == 1]
#' Time.2 = Pooh$Likert[Pooh$Time == 2]
#' library(DescTools)
#' SignTest(x = Time.1, y = Time.2)
#' pairedSampleDominance(x = Time.1, y = Time.2)
#' pairedSampleDominance(Likert ~ Time, data=Pooh)
#' 
#' @importFrom boot boot boot.ci
#' 
#' @export

pairedSampleDominance =
  function(formula=NULL, data=NULL, x=NULL, y=NULL,
           ci=FALSE, conf=0.95, type="perc", R=1000, histogram=FALSE,
           digits=3, na.rm=TRUE, ...)
  {
    
    if(!is.null(formula)){
      z  = eval(parse(text=paste0("data","$",all.vars(formula[[2]])[1])))
      g  = factor(eval(parse(text=paste0("data","$",all.vars(formula[[3]])[1]))))
      x  = z[g==levels(g)[1]]
      y  = z[g==levels(g)[2]]
    }
    
    if(is.null(formula)){
      x = x
      y = y
    }

    if(na.rm){
    Cases = complete.cases(x) & complete.cases(y)
    x = x[Cases]
    y = y[Cases]
    }
    
    Sum    = length(x)
    AdjSum = sum(x[!is.na(x)] != y[!is.na(y)])
    
    if(AdjSum==0){
      Greater = 0
      Less    = 0
      Equal   = 1
    }
    
    if(AdjSum>0){
      Greater = sum(x>y) / Sum
      Less    = sum(x<y) / Sum
      Equal   = sum(x==y) / Sum
    }
    
    es  = Greater - Less
    vda = es / 2 + 0.5
    
    if(ci==TRUE){
      Data = data.frame(x,y)
      Function = function(input, index){
        Input = input[index,]
        ess = (sum(Input$x>Input$y) - sum(Input$x<Input$y)) / length(Input$x)
        return(ess)
      }
      Boot = boot(Data, Function, R=R)
      BCI  = boot.ci(Boot, conf=conf, type=type)
      if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3];}
      if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5];}
      if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5];}
      if(type=="bca") {CI1=BCI$bca[4];      CI2=BCI$bca[5];}
      CI1=signif(CI1, digits=digits)
      CI2=signif(CI2, digits=digits)
      if(is.na(median(x))){CI1=NA; CI2=NA}
      CI3 = CI1 / 2 + 0.5
      CI4 = CI2 / 2 + 0.5
      
      if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray",
                               main="", xlab="dominance")}
      
      Out = data.frame(n=Sum,
                       Less      = signif(Less, digits),
                       Equal     = signif(Equal, digits),
                       Greater   = signif(Greater, digits),
                       Dominance = signif(es, digits),
                       lower.ci  = CI1, 
                       upper.ci  = CI2,
                       VDA      =  signif(vda, digits),
                       lower.vda.ci  = CI3, 
                       upper.vda.ci  = CI4)
    }
    if(ci==FALSE){
      Out = data.frame(n=Sum,
                       Less      = signif(Less, digits),
                       Equal     = signif(Equal, digits),
                       Greater   = signif(Greater, digits),
                       Dominance = signif(es, digits),
                       VDA      =  signif(vda, digits))
    }
    return(Out)
  }