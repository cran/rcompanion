% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fullPTable.r
\name{fullPTable}
\alias{fullPTable}
\title{Convert a lower triangle matrix to a full matrix}
\usage{
fullPTable(PT)
}
\arguments{
\item{PT}{A lower triangle matrix.}
}
\value{
A full matrix.
}
\description{
Converts a lower triangle matrix to a full matrix.
}
\details{
This function is useful to convert a lower triangle matrix
          of p-values from a pairwise test to a full matrix.
          A full matrix can be passed to \code{multcompLetters}
          in the \code{multcompView} package to produce a compact
          letter display.
}
\examples{
### Example with pairwise.wilcox.test
data(BrendonSmall)
BrendonSmall$Instructor = factor(BrendonSmall$Instructor,
                          levels = c('Brendon Small', 'Jason Penopolis',
                                     'Paula Small', 'Melissa Robbins', 
                                     'Coach McGuirk'))
P   = pairwise.wilcox.test(x = BrendonSmall$Score, g = BrendonSmall$Instructor)
PT  = P$p.value
PT
PT1 = fullPTable(PT)
PT1
library(multcompView)
multcompLetters(PT1)

}
\references{
\url{https://rcompanion.org/handbook/F_08.html}
}
\seealso{
\code{\link{cldList}}
}
\author{
Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
}
\concept{cld}
\concept{compact letter display}
\concept{multiple comparisons}
\concept{post-hoc}
