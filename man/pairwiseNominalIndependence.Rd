% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pairwiseNominalIndependence.r
\name{pairwiseNominalIndependence}
\alias{pairwiseNominalIndependence}
\title{Pairwise tests of independence for nominal data}
\usage{
pairwiseNominalIndependence(
  x,
  compare = "row",
  fisher = TRUE,
  gtest = TRUE,
  chisq = TRUE,
  method = "fdr",
  correct = "none",
  yates = FALSE,
  stats = FALSE,
  cramer = FALSE,
  digits = 3,
  ...
)
}
\arguments{
\item{x}{A two-way contingency table. At least one dimension should have
more than two levels.}

\item{compare}{If \code{"row"}, treats the rows as the grouping variable.
If \code{"column"}, treats the columns as the grouping variable.}

\item{fisher}{If \code{"TRUE"}, conducts fisher exact test.}

\item{gtest}{If \code{"TRUE"}, conducts G-test.}

\item{chisq}{If \code{"TRUE"}, conducts Chi-square test of association.}

\item{method}{The method to adjust multiple p-values. 
See \code{stats::p.adjust}.}

\item{correct}{The correction method to pass to \code{DescTools::GTest}.}

\item{yates}{Passed to \code{correct} in \code{stats::chisq.test}.}

\item{stats}{If \code{"TRUE"}, includes the Chi-square value and degrees of 
freedom for Chi-square tests, and the G value.}

\item{cramer}{If \code{"TRUE"}, includes an effect size, Cramer's V in the
output.}

\item{digits}{The number of significant digits in the output.}

\item{...}{Additional arguments, passed to \code{stats::fisher.test}, 
\code{DescTools::GTest}, or \code{stats::chisq.test}.}
}
\value{
A data frame of comparisons, p-values, and adjusted p-values.
}
\description{
Conducts pairwise tests for a 2-dimensional matrix,
             in which at at least one dimension has more than two
             levels, as a post-hoc test. Conducts Fisher exact, Chi-square,
             or G-test.
}
\section{Acknowledgments}{

         My thanks to
         Carole Elliott of Kings Park & Botanic Gardens
         for suggesting the inclusion on the chi-square statistic
         and degrees of freedom in the output.
}

\examples{
### Independence test for a 4 x 2 matrix
data(Anderson)
fisher.test(Anderson)
Anderson = Anderson[(c("Heimlich", "Bloom", "Dougal", "Cobblestone")),]
PT = pairwiseNominalIndependence(Anderson,
                                 fisher = TRUE,
                                 gtest  = FALSE,
                                 chisq  = FALSE,
                                 cramer = TRUE)
PT                                
cldList(comparison = PT$Comparison,
        p.value    = PT$p.adj.Fisher,
        threshold  = 0.05)                             
                                                              
}
\references{
\url{https://rcompanion.org/handbook/H_04.html}
}
\seealso{
\code{\link{pairwiseMcnemar}}, \code{\link{groupwiseCMH}},
          \code{\link{nominalSymmetryTest}}, 
          \code{\link{pairwiseNominalMatrix}}
}
\author{
Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
}
\concept{Fisher's exact test}
\concept{G test}
\concept{chi square test}
\concept{post-hoc}
