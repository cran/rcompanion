% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pairwisePermutationSymmetryMatrix.r
\name{pairwisePermutationSymmetryMatrix}
\alias{pairwisePermutationSymmetryMatrix}
\title{Pairwise two-sample symmetry tests with matrix output}
\usage{
pairwisePermutationSymmetryMatrix(
  formula = NULL,
  data = NULL,
  x = NULL,
  g = NULL,
  b = NULL,
  method = "fdr",
  ...
)
}
\arguments{
\item{formula}{A formula indicating the measurement variable and
the grouping variable. e.g. y ~ group.}

\item{data}{The data frame to use.}

\item{x}{The response variable as a vector.}

\item{g}{The grouping variable as a vector.}

\item{b}{The blocking variable as a vector.}

\item{method}{The p-value adjustment method to use for multiple tests.
See \code{stats::p.adjust}.}

\item{...}{Additional arguments passed to
\code{coin::symmetry_test}.}
}
\value{
A list consisting of:
        A matrix of p-values;
        the p-value adjustment method;
        a matrix of adjusted p-values.
}
\description{
Conducts pairwise two-sample symmetry tests across groups.
}
\details{
The input should include either \code{formula} and \code{data};
         or \code{x}, \code{g}, and \code{b}.
         
         This function is a wrapper for \code{coin::symmetry_test},
         passing pairwise groups to the function. It's critical to read
         and understand the documentation for this function to understand
         its use and options.
}
\note{
The parsing of the formula is simplistic. 
         The first variable on the
         left side is used as the measurement variable.  
         The first variable on the
         right side is used for the grouping variable.
         The second variable on the
         right side is used for the blocking variable.
}
\examples{
data(BobBelcher)

BobBelcher$Instructor = factor( BobBelcher$Instructor, 
                                levels = c("Linda Belcher", "Louise Belcher",
                                           "Tina Belcher", "Bob Belcher",
                                           "Gene Belcher"))

library(coin)

symmetry_test(Likert ~ Instructor | Rater, data= BobBelcher,
              ytrafo   = rank_trafo,
              teststat = "quadratic")

PT = pairwisePermutationSymmetryMatrix(Likert ~ Instructor | Rater,
                                 data     = BobBelcher,
                                 ytrafo   = rank_trafo,
                                 teststat = "quadratic",
                                 method   = "fdr")
PT

PA = PT$Adjusted
library(multcompView)
multcompLetters(PA,
                compare="<",
                threshold=0.05,
                Letters=letters)     

}
\references{
\url{https://rcompanion.org/handbook/K_03.html}
}
\seealso{
\code{\link{pairwisePermutationSymmetry}}
}
\author{
Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
}
\concept{permutation test}
\concept{post-hoc}
