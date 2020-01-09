#' @title Freeman's theta
#'
#' @description Calculates Freeman's theta for a table with one ordinal
#'              variable and one nominal variable; confidence intervals
#'              by bootstrap.
#' 
#' @param x Either a two-way table or a two-way matrix.
#'          Can also be a vector of observations of an ordinal variable.
#' @param g If \code{x} is a vector, \code{g} is the vector of observations for
#'          the grouping, nominal variable.
#' @param group If \code{x} is a table or matrix, \code{group} indicates whether
#'              the \code{"row"} or the \code{"column"} variable is
#'              the nominal, grouping variable.
#' @param verbose If \code{TRUE}, prints statistics for each
#'                   comparison.
#' @param progress If \code{TRUE}, prints a message as each comparison is
#'                 conducted.
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
#' @param reportIncomplete If \code{FALSE} (the default),
#'                         \code{NA} will be reported in cases where there
#'                         are instances of the calculation of the statistic
#'                         failing during the bootstrap procedure.
#' 
#' @details  Freeman's coefficent of differentiation (theta)
#'           is used as a measure of association
#'           for a two-way
#'           table with one ordinal and one nominal variable.
#'           See Freeman (1965).
#'           
#'           Currently, the function makes no provisions for \code{NA}
#'           values in the data.  It is recommended that \code{NA}s be removed
#'           beforehand.
#'           
#'           Because theta is always positive, if \code{type="perc"},
#'           the confidence interval will
#'           never cross zero, and should not
#'           be used for statistical inference.
#'           However, if \code{type="norm"}, the confidence interval
#'           may cross zero.
#'
#'           When theta is close to 0 or very large,
#'           or with small counts in some cells,
#'           the confidence intervals 
#'           determined by this
#'           method may not be reliable, or the procedure may fail.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references Freeman, L.C. 1965. Elementary Applied Statistics for Students
#'             in Behavioral Science. Wiley.
#' 
#'             \url{http://rcompanion.org/handbook/H_11.html}
#' @seealso \code{\link{epsilonSquared}}             
#' @concept correlation Freeman theta ordinal nominal
#' @return A single statistic, Freeman's theta.
#'         Or a small data frame consisting of Freeman's theta,
#'         and the lower and upper confidence limits. 
#'         
#' @examples
#' data(Breakfast)
#' library(coin)
#' chisq_test(Breakfast, scores = list("Breakfast" = c(-2, -1, 0, 1, 2)))
#' freemanTheta(Breakfast)
#' 
#' ### Example from Freeman (1965), Table 10.6
#'Counts = c(1,2,5,2,0,10,5,5,0,0,0,0,2,2,1,0,0,0,2,3)
#'Matrix = matrix(Counts, byrow=TRUE, ncol=5,
#'                dimnames = list(Marital.status=c("Single","Married","Widowed",
#'                                                 "Divorced"),
#'                                Social.adjustment = c("5","4","3","2","1")))
#'Matrix
#'freemanTheta(Matrix)
#' 
#' ### Example after Kruskal Wallis test
#' data(PoohPiglet)
#' kruskal.test(Likert ~ Speaker, data = PoohPiglet)
#' freemanTheta(x = PoohPiglet$Likert, g = PoohPiglet$Speaker)
#' 
#' ### Same data, as table of counts
#' data(PoohPiglet)
#' XT = xtabs( ~ Speaker + Likert , data = PoohPiglet)
#' freemanTheta(XT)
#' 
#' @importFrom boot boot boot.ci
#' 
#' @export

freemanTheta = function (x, g=NULL, group="row", 
                         verbose=FALSE, progress=FALSE,
                         ci=FALSE, conf=0.95, 
                         type="perc",
                         R=1000, histogram=FALSE, digits=3,
                         reportIncomplete=FALSE){
  
  if(is.matrix(x)){x=as.table(x)}
  
  if(is.table(x)){
     Counts = as.data.frame(x)
     Long = Counts[rep(row.names(Counts), Counts$Freq), c(1, 2)]
     rownames(Long) = seq(1:nrow(Long))
     if(group=="row"){
        g=factor(Long[,1])
        x=as.numeric(Long[,2])
        }
     if(group=="column"){
        g=factor(Long[,1])
        x=as.numeric(Long[,2])
     }
  }
  g = factor(g)
  g = droplevels(g)
  k     = length(levels(g))
  Delta = 0
  Tee   = 0
  Count = 0
  Di    = rep(NA, ((k-1)*(k-2)))
  Ti    = rep(NA, ((k-1)*(k-2)))
  for(i in 1:(k-1)){
     for(j in (i+1):k){
       Count = Count + 1
       Y1 = x[g==levels(g)[i]]
       Y2 = x[g==levels(g)[j]]
       n1 = length(Y1)
       n2 = length(Y2)
       tee   = 0
       delta     = 0
       for(l in 1:n1){
         for(m in 1:n2){
           delta = delta + sum(Y1[l] > Y2[m]) - sum(Y1[l] < Y2[m])
         }
       }
       if(progress){cat("Comparison ", Count, " ...\n")}
       Di[Count] = abs(delta)
       Ti[Count] = n1 * n2
       Delta = Delta + Di[Count]
       Tee   = Tee + Ti[Count]
     }
  }

  if(verbose){
  Z = data.frame(Comparison = 1:Count,
                 Di         = Di,
                 Ti         = Ti)
  cat("\n")
  print(Z)
  cat("\n")
  cat("Sum Di = ", Delta,"\n")
  cat("T2     = ", Tee,"\n", "\n")
}
  
  theta = Delta / Tee
  Theta = signif(theta, digits=digits)
  
if(ci==TRUE){
  
    L1     = length(unique(droplevels(factor(g))))
  
  Function = function(input, index){
                      Input = input[index,]
  Input$g = factor(Input$g)
  Input$g = droplevels(Input$g)
  k       = length(levels(Input$g))
  
    NOTEQUAL=0
    if(k != L1){NOTEQUAL=1}
             
    if(NOTEQUAL==1){FLAG=1; return(c(NA,FLAG))}
             
    if(NOTEQUAL==0){
  
       Delta = 0
       Tee   = 0
       Count = 0
       Di    = rep(NA, ((k-1)*(k-2)))
       Ti    = rep(NA, ((k-1)*(k-2)))
        for(i in 1:(k-1)){
          for(j in (i+1):k){
            Count = Count + 1
            Y1 = Input$x[Input$g==levels(Input$g)[i]]
            Y2 = Input$x[Input$g==levels(Input$g)[j]]
            n1 = length(Y1)
            n2 = length(Y2)
            tee   = 0
            delta = 0
              for(l in 1:n1){
                for(m in 1:n2){
                  delta = delta + sum(Y1[l] > Y2[m]) - sum(Y1[l] < Y2[m])
               }}
            Di[Count] = abs(delta)
            Ti[Count] = n1 * n2
            Delta = Delta + Di[Count]
            Tee = Tee + Ti[Count]
         }
        }
    FLAG = 0   
    return(c((Delta / Tee),FLAG))}
  }
  Data = data.frame(x,g)
  Boot = boot(Data, Function, R=R)
  BCI  = boot.ci(Boot, conf=conf, type=type)
  if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3];}
  if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5];}
  if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5];}
  if(type=="bca")  {CI1=BCI$bca[4];     CI2=BCI$bca[5];}  
  
  if(sum(Boot$t[,2])>0 & reportIncomplete==FALSE) {CI1=NA; CI2=NA}
  
  CI1=signif(CI1, digits=digits)
  CI2=signif(CI2, digits=digits)
  
  if(histogram==TRUE){hist(Boot$t[,1], 
                      col = "darkgray", xlab="theta", main="")}
  
}
  
if(ci==FALSE){names(Theta)="Freeman.theta"; return(Theta)}
if(ci==TRUE){names(Theta) = ""
             return(data.frame(Freeman.theta=Theta, lower.ci=CI1, upper.ci=CI2))}
}
