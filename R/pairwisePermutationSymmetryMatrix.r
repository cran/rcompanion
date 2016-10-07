#' @title Pairwise two-sample permutation symmetry 
#'        tests with matrix output
#'
#' @description Conducts pairwise two-sample permutation tests
#'              for symmetry across groups.
#' 
#' @param x      The response variable as a vector.
#' @param g      The grouping variable as a vector.
#' @param b      The blocking variable as a vector.
#' @param method The p-value adjustment method to use for multiple tests.
#'               See \code{\link{p.adjust}}.
#' @param ...    Additional arguments passed to
#'               \code{\link{symmetry_test}}.               
#'             
#' @details Permutation tests are non-parametric tests 
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
#' @examples
#' data(BobBelcher)
#' BobBelcher = BobBelcher[order(factor(BobBelcher$Instructor, 
#'                         levels=c("Linda Belcher", "Louise Belcher",
#'                                  "Tina Belcher", "Bob Belcher",
#'                                  "Gene Belcher"))),]
#' BobBelcher$Likert.f = factor(BobBelcher$Likert, ordered=TRUE)
#' PT = pairwisePermutationSymmetryMatrix(x      = BobBelcher$Likert.f,
#'                                        g      = BobBelcher$Instructor,
#'                                        b      = BobBelcher$Rater,
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
  function(x, g, b, method = "fdr", ...)
  {
  n = length(unique(g))
  N = n*n
  d = data.frame(x = x, g = g, b = b)
  Y = matrix(rep(NA_real_, N),ncol=n)
  rownames(Y)=unique(g)
  colnames(Y)=unique(g)
  Z = matrix(rep(NA_real_, N),ncol=n)
  rownames(Z)=unique(g)
  colnames(Z)=unique(g)
  k=0
  for(i in 1:(n-1)){
     for(j in (i+1):n){
     k=k+1
     Datax = subset(d, g==unique(g)[i])
     Datay = subset(d, g==unique(g)[j])
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