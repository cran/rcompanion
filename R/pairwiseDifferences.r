#' @title Pairwise differences for unreplicated CBD 
#'
#' @description Calculates the differences in the response variable
#'              for each pair of levels of a grouping variable 
#'              in an unreplicated complete block design.
#'              
#' @param formula A formula indicating the measurement variable and
#'                the grouping variable. e.g. y ~ group.
#' @param data   The data frame to use.
#' @param x The vector of the response variable.
#' @param g The vector of the grouping variable.
#' @param plotit If \code{TRUE}, then produces bar plots of the differences.
#' @param factorize If \code{TRUE}, then adds a column to the output 
#'                  data frame consisting of the differences as
#'                  a factor variable. This output is added automatically if
#'                  \code{plotit = TRUE}.          
#'             
#' @details The main use of the function is to check the shape of the
#'          distribution of differences in responses for paired t-test,
#'          paired rank-sum test, Friedman test, or Quade test.
#' 
#'           The function assumes that the data frame is already ordered by
#'           the blocking variable, so that the first observation of Group 1
#'           is paired with the first observation of Group 2, and so on.
#'           
#'           The function assumes that the data are in complete block design.
#'           That is,
#'           for any level of the grouping variable in Group 1 there exists
#'           one paired value in Group 2, and so on.
#'           
#'           The input should include either \code{formula} and \code{data};
#'           or \code{x}, and \code{g}.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_10.html}
#' @concept paired t-test rank-sum friedman quade histogram
#' @return A data frame of the paired groups, the differences in
#'         their response variables,
#'         and optionally the differences expressed as a factor variable.
#'         If \code{plotit = TRUE}, then also produce one or more plots.
#'
#' @note    The parsing of the formula is simplistic. 
#'          The first variable on the
#'          left side is used as the measurement variable.  
#'          The first variable on the
#'          right side is used for the grouping variable.
#'          
#' @examples
#' ### Two-sample paired data example with bar plot
#' data(Pooh)
#' Pooh.diff = pairwiseDifferences(Likert ~ Time, 
#'                                 data=Pooh, 
#'                                 plotit = TRUE)
#' 
#' ### Unreplicated complete block design example with bar plots
#' data(BobBelcher)
#' Bob.diff = pairwiseDifferences(Likert ~ Instructor,
#'                                data=BobBelcher, 
#'                                factorize=TRUE)
#' library(lattice)
#' histogram(~ Difference.f | Comparison,
#'          data=Bob.diff,
#'          type = "count",
#'          layout=c(2,5))
#' 
#' @importFrom graphics barplot
#' @importFrom lattice  histogram
#' 
#' @export

pairwiseDifferences = 
  function(formula=NULL, data=NULL,
           x=NULL, g=NULL, plotit=FALSE, factorize=FALSE)
  {
  if(!is.null(formula)){
    x  = eval(parse(text=paste0("data","$",all.vars(formula[[2]])[1])))
    g  = eval(parse(text=paste0("data","$",all.vars(formula[[3]])[1])))
    }
  if(!is.factor(g)){g=factor(g)}
  n = length(levels(g))
  N = n*(n-1)/2
  Flength = rep(NA, n)
  i=0
  for(i in (i+1):n){
     Flength[i] = length(g[g==levels(g)[i]])
     }
  m = min(Flength)
  M = N * m
  d = data.frame(x = x, g = g)
  Z = data.frame(Comparison=rep("A", M),
                 Difference=rep(as.numeric(NA), M),
                 Difference.f=rep("A", M),
                 stringsAsFactors=FALSE)
  i=0; j=0; k=0; l=0              
  for(i in 1:(n-1)){
     for(j in (i+1):n){
     Namea = as.character(levels(g)[i])
     Nameb = as.character(levels(g)[j])
     Datax = subset(d, g==levels(g)[i])
     Datay = subset(d, g==levels(g)[j])
     for(k in 1:m){
     l=l+1
     a = Datax$x[k]
     b = Datay$x[k]
     Diff = as.numeric(a-b)
     Z[l,] = c( paste0(Namea, " - ", Nameb), 
                as.numeric(Diff),
                "NA")
       }
      }
     }
  Z$Difference = as.numeric(Z$Difference)
  if(factorize|plotit){
     a=min(Z$Difference)
     b=max(Z$Difference)
     Frange = as.character(rep(a:b))
     Z$Difference.f = factor(Z$Difference,
                      ordered = TRUE,
                      levels=Frange)
     }
  if(plotit){
     for(i in 1:(N)){
     Y = subset(Z, Z$Comparison==unique(Z$Comparison)[i])
     X = xtabs(~ Y$Difference.f)
     Name=Y$Comparison[1]
        barplot(X,   
                col="dark gray", 
                xlab=bquote(paste(.(Name), " | Difference")),
                ylab="Frequency")
     }
    }
  Z
  }