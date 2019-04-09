#' @title Cramer's V for chi-square goodness-of-fit tests
#'
#' @description Calculates Cramer's V for a vector of counts and expected 
#'              counts; confidence intervals by bootstrap.
#' 
#' @param x A vector of observed counts.
#' @param p A vector of expected or default probabilities.
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
#' @param ...    Additional arguments passed to \code{chisq.test}. 
#' 
#' @details This modification of Cramer's V
#'          could be used to indicate an effect size
#'          in cases where a chi-square goodness-of-fit test might be used.
#'          It indicates the degree of deviation of observed counts
#'          from the expected probabilities.
#'           
#'           In the case of equally-distributed expected frequencies,
#'           Cramer's V will be equal to 1 when all counts are in one category,
#'           and it will be equal to 0 when the counts are equally distributed
#'           across categories.
#'           This does not hold if the expected frequencies are not
#'           equally-distributed. 
#'           
#'           When V is close to 0 or 1,
#'           or with small counts, 
#'           the confidence intervals 
#'           determined by this
#'           method may not be reliable, or the procedure may fail.
#' 
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/H_03.html}
#' @concept correlation effect size cramer V
#' @seealso \code{\link{cramerV}}
#' @return A single statistic, Cramer's V.  
#'         Or a small data frame consisting of Cramer's V,
#'         and the lower and upper confidence limits.
#'         
#' @examples
#' ### Equal probabilities example
#' ### From http://rcompanion.org/handbook/H_03.html
#' nail.color = c("Red", "None", "White", "Green", "Purple", "Blue")
#' observed   = c( 19,    3,      1,       1,       2,        2    )
#' expected   = c( 1/6,   1/6,    1/6,     1/6,     1/6,      1/6  )
#' chisq.test(x = observed, p = expected)
#' cramerVFit(x = observed, p = expected)
#' 
#' ### Unequal probabilities example
#' ### From http://rcompanion.org/handbook/H_03.html
#' race = c("White", "Black", "American Indian", "Asian", "Pacific Islander",
#'           "Two or more races")
#' observed = c(20, 9, 9, 1, 1, 1)
#' expected = c(0.775, 0.132, 0.012, 0.054, 0.002, 0.025)
#' chisq.test(x = observed, p = expected)
#' cramerVFit(x = observed, p = expected)
#' 
#' ### Examples of perfect and zero fits
#' cramerVFit(c(100, 0, 0, 0, 0))
#' cramerVFit(c(10, 10, 10, 10, 10))
#' 
#' @importFrom stats chisq.test
#' @importFrom boot  boot boot.ci
#' 
#' @export

cramerVFit = function(x, p=rep(1/length(x), length(x)), 
                      ci=FALSE, conf=0.95, type="perc",
                      R=1000, histogram=FALSE, digits=4, ...) {
  CV=NULL
  N = sum(x)
  Chi.sq = suppressWarnings(chisq.test(x=x, p=p, ...)$statistic)
  K   = length(x)
  CV =  sqrt(Chi.sq/N/(K-1))
  
  CV = signif(as.numeric(CV), digits=digits)
  
if(ci==TRUE){
    Counts = as.data.frame(x)
    Long=data.frame(x = rep(row.names(Counts), Counts$x))
    rownames(Long) = seq(1:nrow(Long))
    
  Function = function(input, index){
             Input = input[index,]
             
  Obs = as.vector(table(Input))
  CV=NULL
  N = sum(Obs)
  Chi.sq = suppressWarnings(chisq.test(x=Obs, p=p, ...)$statistic)
  K   = length(Obs)
  CV =  sqrt(Chi.sq/N/(K-1))
  
  CV = signif(as.numeric(CV), digits=digits)
  
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
  
  if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray",
                      main="", xlab="V")}
  
}
 if(ci==FALSE){names(CV)="Cramer V"; return(CV)}
 if(ci==TRUE){return(data.frame(r=CV, lower.ci=CI1, upper.ci=CI2))}  
}
