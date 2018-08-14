#' @title Freeman's theta
#'
#' @description Calculates Freeman's theta for a table with one ordinal
#'              variable and one nominal variable.
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
#' @param digits The number of significant digits in the output.
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
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references Freeman, L.C. 1965. Elementary Applied Statistics for Students
#'             in Behavioral Science. Wiley.
#' 
#'             \url{http://rcompanion.org/handbook/H_11.html}
#' @seealso \code{\link{epsilonSquared}}             
#' @concept correlation Freeman theta ordinal nominal
#' @return A single statistic, Freeman's theta
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
#' @export

freemanTheta = function (x, g=NULL, group="row", 
                         verbose=FALSE, progress=FALSE,
                         digits=3){
  
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
  g=droplevels(g)
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
  names(Theta)="Freeman.theta"

  return(Theta)
}
