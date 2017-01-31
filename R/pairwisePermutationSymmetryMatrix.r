#' @title Pairwise two-sample permutation symmetry 
#'        tests with matrix output
#'
#' @description Conducts pairwise two-sample permutation tests
#'              for symmetry across groups.
#'              
#' @param formula A formula indicating the measurement variable and
#'                the grouping variable. e.g. y ~ group.
#' @param data   The data frame to use.
#' @param x      The response variable as a vector.
#' @param g      The grouping variable as a vector.
#' @param b      The blocking variable as a vector.
#' @param method The p-value adjustment method to use for multiple tests.
#'               See \code{\link{p.adjust}}.
#' @param ...    Additional arguments passed to
#'               \code{\link{symmetry_test}}.               
#'             
#' @details The input should include either \code{formula} and \code{data};
#'          or \code{x}, \code{g}, and \code{b}.
#'          
#'          Permutation tests are non-parametric tests 
#'          that do not assume normally-distributed errors.
#'          See \url{http://rcompanion.org/rcompanion/d_06a.html} for
#'          futher discussion of this test.
#' 
#'          The \code{pairwisePermutationSymmetryMatrix} function
#'          can be used as a post-hoc method following an omnibus 
#'          permutation test analogous to a paired one-way analysis
#'          of variance.
#'          The matrix output can be converted to a compact letter display.                                                                                   
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/K_03.html}
#' @seealso \code{\link{pairwisePermutationSymmetry}}
#' @concept permutation nonparametric post-hoc one-way cld symmetry
#' @return A list consisting of:
#'         A matrix of p-values;
#'         the p-value adjustment method;
#'         a matrix of adjusted p-values. 
#'      
#' @note    The parsing of the formula is simplistic. 
#'          The first variable on the
#'          left side is used as the measurement variable.  
#'          The first variable on the
#'          right side is used for the grouping variable.
#'          The second variable on the
#'          right side is used for the blocking variable.   
#'         
#' @examples
#' data(BobBelcher)
#' BobBelcher$Instructor = factor( BobBelcher$Instructor, 
#'                   levels = c("Linda Belcher", "Louise Belcher",
#'                              "Tina Belcher", "Bob Belcher",
#'                              "Gene Belcher"))
#' BobBelcher$Likert.f = factor(BobBelcher$Likert, ordered=TRUE)
#' PT = pairwisePermutationSymmetryMatrix(Likert.f ~ Instructor | Rater,
#'                                        data   = BobBelcher,
#'                                        method = "fdr")$Adjusted
#' PT
#' library(multcompView)
#' multcompLetters(PT,
#'                 compare="<",
#'                 threshold=0.05,
#'                 Letters=letters)   
#' 
#' @importFrom stats p.adjust
#' @importFrom coin symmetry_test statistic
#' 
#' @export  

pairwisePermutationSymmetryMatrix = 
  function(formula=NULL, data=NULL,
           x=NULL, g=NULL, b=NULL, 
           method = "fdr", ...)
  {
  if(!is.null(formula)){
    x  = eval(parse(text=paste0("data","$",all.vars(formula[[2]])[1])))
    g  = eval(parse(text=paste0("data","$",all.vars(formula[[3]])[1])))
    b  = eval(parse(text=paste0("data","$",all.vars(formula[[3]])[2])))
    }
  if(!is.factor(g)){g=factor(g)}
  if(!is.factor(b)){g=factor(b)}
  n = length(levels(g))
  N = n*n
  d = data.frame(x = x, g = g, b = b)
  Y = matrix(rep(NA_real_, N),ncol=n)
  rownames(Y)=levels(g)
  colnames(Y)=levels(g)
  Z = matrix(rep(NA_real_, N),ncol=n)
  rownames(Z)=levels(g)
  colnames(Z)=levels(g)
  k=0
  for(i in 1:(n-1)){
     for(j in (i+1):n){
     k=k+1
     Datax = subset(d, g==levels(g)[i])
     Datay = subset(d, g==levels(g)[j])
     Dataz = rbind(Datax, Datay)
     Dataz$g2 = factor(Dataz$g)
     z = symmetry_test(x ~ g2|b, data=Dataz, ...)                  
   Y[i,j] = signif(pvalue(z), digits = 4)
   } 
   }
Z[upper.tri(Z)] = 
      signif(p.adjust(Y[upper.tri(Y)], method=method), digits=4)
Z = t(Z)
Z[upper.tri(Z)] = 
      signif(p.adjust(Y[upper.tri(Y)], method=method), digits=4)
diag(Z) = signif(1.00, digits = 4)
W = method
V = list(Y, W, Z)
names(V) = c("Unadjusted",
             "Method",
             "Adjusted")
return(V)   
} 