#' @title Cramer's V (phi)
#'
#' @description Calculates Cramer's V for a table of nominal variables;
#'              confidence intervals by bootstrap.
#' 
#' @param x Either a two-way table or a two-way matrix.
#'          Can also be a vector of observations for one dimension
#'          of a two-way table. 
#' @param y If \code{x} is a vector, \code{y} is the vector of observations for
#'          the second dimension of a two-way table.
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
#' @param bias.correct If \code{TRUE}, a bias correction is applied.
#' @param ...    Additional arguments passed to \code{chisq.test}. 
#' 
#' @details  Cramer's V is used as a measure of association
#'           between two nominal variables, or as an effect size
#'           for a chi-square test of association.  For a 2 x 2 table,
#'           the absolute value of the phi statistic is the same as
#'           Cramer's V.
#'           
#'           Because V is always positive, if \code{type="perc"},
#'           the confidence interval will
#'           never cross zero. In this case, 
#'           the confidence interval range should not
#'           be used for statistical inference.
#'           However, if \code{type="norm"}, the confidence interval
#'           may cross zero.  
#'           
#'           When V is close to 0 or very large,
#'           or with small counts, 
#'           the confidence intervals 
#'           determined by this
#'           method may not be reliable, or the procedure may fail.
#'                      
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/H_10.html}
#' @concept correlation phi cramer V
#' @seealso \code{\link{cohenW}}
#' @return A single statistic, Cramer's V.  
#'         Or a small data frame consisting of Cramer's V,
#'         and the lower and upper confidence limits.
#'         
#' @examples
#' ### Example with table
#' data(Anderson)
#' fisher.test(Anderson)
#' cramerV(Anderson)
#'
#' ### Example with two vectors
#' Species = c(rep("Species1", 16), rep("Species2", 16))
#' Color   = c(rep(c("blue", "blue", "blue", "green"),4),
#'             rep(c("green", "green", "green", "blue"),4))
#' fisher.test(Species, Color)
#' cramerV(Species, Color)
#' 
#' @importFrom stats chisq.test
#' @importFrom boot  boot boot.ci 
#' @export

cramerV = function(x, y=NULL, 
                   ci=FALSE, conf=0.95, type="perc",
                   R=1000, histogram=FALSE, 
                   digits=4, bias.correct=FALSE, ...) {
  
  CV=NULL
  
  if(is.factor(x)){x=as.vector(x)}
  if(is.factor(y)){x=as.vector(y)}
  if(is.vector(x) & is.vector(y)){
  N      = length(x)
  Chi.sq = suppressWarnings(chisq.test(x, y, correct=FALSE, ...)$statistic)
  Phi    = Chi.sq / N
  Row    = length(unique(x))
  C      = length(unique(y))
  CV     =  sqrt(Phi / min(Row-1, C-1))
  }
  
  if(is.matrix(x)){x=as.table(x)}
  if(is.table(x)){
  N = sum(x)
  Chi.sq = suppressWarnings(chisq.test(x, correct=FALSE, ...)$statistic)
  Phi = Chi.sq / N
  Row = nrow(x)
  C   = ncol(x)
  CV =  sqrt(Phi / min(Row-1, C-1))}
  
  if(bias.correct==TRUE){Phi = max(0, Phi-((Row-1)*(C-1)/(N-1)))
                        CC  = C-((C-1)^2/(N-1))
                        RR  = Row-((Row-1)^2/(N-1))
                        CV  = sqrt(Phi / min(RR-1, CC-1))} 
 
  CV = signif(as.numeric(CV), digits=digits)
  
if(ci==TRUE){
  if(is.matrix(x)){x=as.table(x)}
  if(is.table(x)){
    Counts = as.data.frame(x)
    Long = Counts[rep(row.names(Counts), Counts$Freq), c(1, 2)]
    rownames(Long) = seq(1:nrow(Long))
    }
  if(is.vector(x) & is.vector(y)){  
    Long = data.frame(x=x, y=y)
  } 
  Function = function(input, index){
             Input = input[index,]
  N      = length(Input[,1])
  Chi.sq = suppressWarnings(chisq.test(Input[,1], Input[,2], correct=FALSE, ...)$statistic)
  Phi    =  Chi.sq / N
  Row      = length(unique(Input[,1]))
  C      = length(unique(Input[,2]))
  CV     =  sqrt(Phi / min(Row-1, C-1))
  
  if(bias.correct==TRUE){Phi = max(0, Phi-((Row-1)*(C-1)/(N-1)))
                        CC  = C-((C-1)^2/(N-1))
                        RR  = Row-((Row-1)^2/(N-1))
                        CV  = sqrt(Phi / min(RR-1, CC-1))}
  return(CV)
  }

  Boot = boot(Long, Function, R=R)
  BCI  = boot.ci(Boot, conf=conf, type=type)
  if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3];}
  if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5];}
  if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5];}
  if(type=="bca")  {CI1=BCI$bca[4];     CI2=BCI$bca[5];}  
  
  CI1=signif(CI1, digits=digits)
  CI2=signif(CI2, digits=digits)
  
  if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray", xlab="V", main="")}

}
 if(ci==FALSE){names(CV)="Cramer V"; return(CV)}
 if(ci==TRUE){return(data.frame(Cramer.V=CV, lower.ci=CI1, upper.ci=CI2))}  
}
