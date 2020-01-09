#' @title Cohen's w (omega)
#'
#' @description Calculates Cohen's w for a table of nominal variables.
#' 
#' @param x Either a two-way table or a two-way matrix.
#'          Can also be a vector of observations for one dimension
#'          of a two-way table. 
#' @param y If \code{x} is a vector, \code{y} is the vector of observations for
#'          the second dimension of a two-way table.
#' @param p If \code{x} is a vector of observed counts, \code{p} can be given as
#'          a vector of theoretical probabilties,
#'          as in a chi-square goodness of fit test.
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
#'                         In the case of the goodness-of-fit
#'                         scenario, setting this to \code{TRUE}
#'                         will have no effect.          
#' @param ...    Additional arguments passed to \code{chisq.test}.
#' 
#' @details  Cohen's w is used as a measure of association
#'           between two nominal variables, or as an effect size
#'           for a chi-square test of association.  For a 2 x 2 table,
#'           the absolute value of the phi statistic is the same as
#'           Cohen's w.  
#'           The value of Cohen's w is not bound by 1 on the upper end.
#'
#'           Cohen's w is "naturally nondirectional". That is,
#'           the value will always be zero or positive.
#'           Because of this, if \code{type="perc"},
#'           the confidence interval will
#'           never cross zero.
#'           The confidence interval range should not
#'           be used for statistical inference.
#'           However, if \code{type="norm"}, the confidence interval
#'           may cross zero.  
#'           
#'           When w is close to 0 or very large,
#'           or with small counts, 
#'           the confidence intervals 
#'           determined by this
#'           method may not be reliable, or the procedure may fail.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/H_10.html}
#' @seealso \code{\link{cramerV}} \code{\link{cramerVFit}}
#' @concept correlation phi cohen w omega
#' @return A single statistic, Cohen's w.
#'         Or a small data frame consisting of Cohen's w,
#'         and the lower and upper confidence limits.
#'         
#' @examples
#' ### Example with table
#' data(Anderson)
#' fisher.test(Anderson)
#' cohenW(Anderson)
#' 
#' ### Example for goodness-of-fit
#' ### Bird foraging example, Handbook of Biological Statistics
#' observed = c(70,   79,   3,    4)
#' expected = c(0.54, 0.40, 0.05, 0.01)
#' chisq.test(observed, p = expected)
#' cohenW(observed, p = expected)
#'
#' ### Example with two vectors
#' Species = c(rep("Species1", 16), rep("Species2", 16))
#' Color   = c(rep(c("blue", "blue", "blue", "green"),4),
#'             rep(c("green", "green", "green", "blue"),4))
#' fisher.test(Species, Color)
#' cohenW(Species, Color)
#' 
#' @importFrom stats chisq.test
#' 
#' @export

cohenW = function(x, y=NULL, p=NULL,
                  ci=FALSE, conf=0.95, type="perc",
                  R=1000, histogram=FALSE, 
                  digits=4, reportIncomplete=FALSE, ...) {
  CW=NULL
  if(is.factor(x)){x=as.vector(x)}
  if(is.factor(y)){y=as.vector(y)}
  if(is.vector(x) & is.vector(y)){
  Chi.sq = suppressWarnings(chisq.test(x, y, correct=FALSE, ...))
  }
  
  if(is.vector(x) & !is.null(p)){
  Chi.sq = suppressWarnings(chisq.test(x=x, p=p, correct=FALSE, ...))
  }
  
 if(is.matrix(x)){x=as.table(x)}
  
 if(is.table(x)){
  Chi.sq = suppressWarnings(chisq.test(x, correct=FALSE, ...))
  }
  
  Sum      = sum(Chi.sq$observed)
  Expected = Chi.sq$expected/Sum
  Observed = Chi.sq$observed/Sum

  CW       = sqrt(sum((Observed-Expected)^2/Expected))
  
  CW = signif(as.numeric(CW), digits=digits)
  if(ci==FALSE){names(CW) = "Cohen w"; return(CW)}
  
  if(is.nan(CW) & ci==TRUE){
    return(data.frame(Cohen.w=CW, lower.ci=NA, upper.ci=NA))} 
  
  if(ci==TRUE){
  if(is.matrix(x)){x=as.table(x)}
  if(is.table(x)){
    Type = 1
    Counts = as.data.frame(x)
    Long = Counts[rep(row.names(Counts), Counts$Freq), c(1, 2)]
    rownames(Long) = seq(1:nrow(Long))
    }
  if(is.vector(x) & is.vector(y)){
    Type = 1
    Long = data.frame(x=x, y=y)
  }
    if(is.vector(x) & !is.null(p)){
      Type = 2
      Counts = data.frame(Count = x, Cat = letters[1:length(x)])
      Long = data.frame(Cat = Counts[rep(row.names(Counts), Counts$Count),
              c("Cat")])
      rownames(Long) = seq(1:nrow(Long))
    }
    
  if(Type==1){  
    L1     = length(unique(droplevels(Long[,1])))
    L2     = length(unique(droplevels(Long[,2])))
  }
    
    if(Type==2){
    L1     = length(unique(droplevels(Long$Cat)))
  }  
    
  Function = function(input, index){
    Input = input[index,]
             
    NOTEQUAL=0
    if(Type==1){
      if(length(unique(droplevels(Input[,1]))) != L1 |
                length(unique(droplevels(Input[,2]))) != L2){NOTEQUAL=1}}
    
    if(Type==2){
      if(length(unique(droplevels(Input))) != L1){NOTEQUAL=1}}
             
    if(NOTEQUAL==1){FLAG=1; return(c(NA,FLAG))}
             
    if(NOTEQUAL==0){
             
    if(Type==1){
       Chi.sq = suppressWarnings(chisq.test(Input[,1], Input[,2], 
         correct=FALSE, ...))
       }
  
    if(Type==2){
       Chi.sq = suppressWarnings(chisq.test(x=table(Input), p=p, 
         correct=FALSE, ...))
       }
  
    Sum      = sum(Chi.sq$observed)
    Expected = Chi.sq$expected/Sum
    Observed = Chi.sq$observed/Sum

    CW       = sqrt(sum((Expected-Observed)^2/Expected))
    FLAG     = 0
    return(c(CW,FLAG))}
    }

  Boot = boot(Long, Function, R=R)
  BCI  = boot.ci(Boot, conf=conf, type=type)
  if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3]}
  if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5]}
  if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5]}
  if(type=="bca")  {CI1=BCI$bca[4];     CI2=BCI$bca[5]}
  
  if(Type==1 & sum(Boot$t[,2])>0 & reportIncomplete==FALSE) {CI1=NA; CI2=NA}
  if(Type==2 & sum(Boot$t[,2])>0) {CI1=NA; CI2=NA}
  
  CI1=signif(CI1, digits=digits)
  CI2=signif(CI2, digits=digits)
  
  if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray", xlab="w", main="")}

}
 if(ci==TRUE){return(data.frame(Cohen.w=CW, lower.ci=CI1, upper.ci=CI2))}  
}

