#' @title phi
#'
#' @description Calculates phi for a 2 x 2 table of nominal variables;
#'              confidence intervals by bootstrap.
#' 
#' @param x Either a 2 x 2 table or a 2 x 2 matrix.
#'          Can also be a vector of observations for one dimension
#'          of a 2 x 2 table. 
#' @param y If \code{x} is a vector, \code{y} is the vector of observations for
#'          the second dimension of a 2 x2  table.
#' @param ci If \code{TRUE}, returns confidence intervals by bootstrap.
#'           May be slow.
#' @param conf The level for the confidence interval.
#' @param type The type of confidence interval to use.
#'             Can be any of "\code{norm}", "\code{basic}", 
#'                           "\code{perc}", or "\code{bca}".
#'             Passed to \code{boot.ci}.
#' @param R The number of replications to use for bootstrap.
#' @param histogram If \code{TRUE}, produces a histogram of bootstrapped values.
#' @param verbose If \code{TRUE}, prints the table of counts.
#' @param digits The number of significant digits in the output.
#' @param reportIncomplete If \code{FALSE} (the default),
#'                         \code{NA} will be reported in cases where there
#'                         are instances of the calculation of the statistic
#'                         failing during the bootstrap procedure.             
#' @param ...    Additional arguments. (Ignored.) 
#' 
#' @details  phi is used as a measure of association
#'           between two binomial variables, or as an effect size
#'           for a chi-square test of association for a 2 x 2 table.
#'           The absolute value of the phi statistic is the same as
#'           Cramer's V for a 2 x 2 table.
#'           
#'           Unlike Cramer's V, phi can be positive or negative (or zero), and
#'           ranges from -1 to 1.
#'           
#'           When phi is close to its extremes,
#'           or with small counts, 
#'           the confidence intervals 
#'           determined by this
#'           method may not be reliable, or the procedure may fail.
#'                      
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/H_10.html}
#' @concept correlation phi cramer V
#' @seealso \code{\link{cramerV}}
#' @return A single statistic, phi.  
#'         Or a small data frame consisting of phi,
#'         and the lower and upper confidence limits.
#'         
#' @examples
#' ### Example with table
#' Matrix = matrix(c(13, 26, 26, 13), ncol=2)
#' phi(Matrix)
#'
#' ### Example with two vectors
#' Species = c(rep("Species1", 16), rep("Species2", 16))
#' Color   = c(rep(c("blue", "blue", "blue", "green"),4),
#'             rep(c("green", "green", "green", "blue"),4))
#' phi(Species, Color)
#' 
#' @importFrom boot  boot boot.ci 
#' @export

phi = function(x, y=NULL, 
                   ci=FALSE, conf=0.95, type="perc",
                   R=1000, histogram=FALSE, verbose=FALSE, 
                   digits=3,
                   reportIncomplete=FALSE, ...) {
  
  PHI=NULL
  
  if(is.factor(x)){x=as.vector(x)}
  if(is.factor(y)){x=as.vector(y)}
  if(is.vector(x) & is.vector(y)){
  if((length(unique(x)) != 2) || (length(unique(y)) != 2))
     {stop("phi is applicable only for 2 binomial variables")}
     Tab = xtabs(~ x + y)
  }
  
  if(is.matrix(x)){Tab=as.table(x)}
  if(is.table(x)){Tab = x}
    if((nrow(Tab) != 2) || (ncol(Tab) != 2))
    {stop("phi is applicable only for a 2 x 2 table")}
  if(verbose){print(Tab) ;cat("\n")}
  Tab2 = Tab / sum(Tab)
  a=Tab2[1,1]; b=Tab2[1,2]; c=Tab2[2,1]; d=Tab2[2,2]
  PHI = (a- (a+b)*(a+c))/sqrt((a+b)*(c+d)*(a+c)*(b+d) )
  Phi= signif(as.numeric(PHI), digits=digits)
  
  if(is.nan(Phi) & ci==TRUE){
    return(data.frame(phi=Phi, lower.ci=NA, upper.ci=NA))} 
  
if(ci==TRUE){
    Counts = as.data.frame(Tab)
    Long = Counts[rep(row.names(Counts), Counts$Freq), c(1, 2)]
    rownames(Long) = seq(1:nrow(Long))
    
    L1     = length(unique(droplevels(Long[,1])))
    L2     = length(unique(droplevels(Long[,2])))
    
  Function = function(input, index){
             Input = input[index,]
             
             NOTEQUAL=0
             if(length(unique(droplevels(Input[,1]))) != L1 |
                length(unique(droplevels(Input[,2]))) != L2){NOTEQUAL=1}
             
             if(NOTEQUAL==1){FLAG=1; return(c(NA,FLAG))}
             
             if(NOTEQUAL==0){
                Tab = xtabs(~ Input[,1] + Input[,2])
                Tab2 = Tab / sum(Tab)
                a=Tab2[1,1]; b=Tab2[1,2]; c=Tab2[2,1]; d=Tab2[2,2]
                PHI = (a- (a+b)*(a+c))/sqrt((a+b)*(c+d)*(a+c)*(b+d))
                FLAG = 0
                return(c(PHI,FLAG))}
             }

  Boot = boot(Long, Function, R=R)
  BCI  = boot.ci(Boot, conf=conf, type=type)
  if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3]}
  if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5]}
  if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5]}
  if(type=="bca")  {CI1=BCI$bca[4];     CI2=BCI$bca[5]}  
  
  if(sum(Boot$t[,2])>0 & reportIncomplete==FALSE) {CI1=NA; CI2=NA}
  
  CI1=signif(CI1, digits=digits)
  CI2=signif(CI2, digits=digits)
  
  if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray", xlab="phi", main="")}

}
 if(ci==FALSE){names(Phi)="phi"; return(Phi)}
 if(ci==TRUE){return(data.frame(phi=Phi, lower.ci=CI1, upper.ci=CI2))}  
}
