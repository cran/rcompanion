% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wilcoxonR.r
\name{wilcoxonR}
\alias{wilcoxonR}
\title{r effect size for Wilcoxon two-sample rank-sum test}
\usage{
wilcoxonR(x, g = NULL, group = "row", digits = 3, ...)
}
\arguments{
\item{x}{Either a two-way table or a two-way matrix.
Can also be a vector of observations of an ordinal variable.}

\item{g}{If \code{x} is a vector, \code{g} is the vector of observations for
the grouping, nominal variable.
Only the first two levels of the nominal variable are used.}

\item{group}{If \code{x} is a table or matrix, \code{group} indicates whether
the \code{"row"} or the \code{"column"} variable is
the nominal, grouping variable.}

\item{digits}{The number of significant digits in the output.}

\item{...}{Additional arguments passed to the \code{wilcox_test} function.}
}
\value{
A single statistic, r
}
\description{
Calculates r effect size
             for Mann-Whitney, two-sample rank-sum test,
             or a table with an ordinal variable and a
             nominal variable with two levels.
}
\details{
A Z value is extracted from the \code{wilcox_test} function in the
          coin package.  r  is calculated as Z divided by 
          square root of the total observations.
 
          Currently, the function makes no provisions for \code{NA}
          values in the data.  It is recommended that \code{NA}s be removed
          beforehand.
}
\examples{
data(Breakfast)
Table = Breakfast[1:2,]
library(coin)
chisq_test(Table, scores = list("Breakfast" = c(-2, -1, 0, 1, 2)))
wilcoxonR(Table)

data(PoohPiglet)
Data = PoohPiglet[PoohPiglet$Speaker \%in\% c("Pooh", "Piglet"),]
wilcox.test(Likert ~ Speaker, data = Data)
wilcoxonR(x = Data$Likert, g = Data$Speaker)

}
\references{
\url{http://rcompanion.org/handbook/F_04.html}
}
\seealso{
\code{\link{freemanTheta}}
}
\author{
Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
}
\concept{
correlation effect size ordinal nominal
}