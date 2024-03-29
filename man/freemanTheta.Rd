% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/freemanTheta.r
\name{freemanTheta}
\alias{freemanTheta}
\title{Freeman's theta}
\usage{
freemanTheta(
  x,
  g = NULL,
  group = "row",
  verbose = FALSE,
  progress = FALSE,
  ci = FALSE,
  conf = 0.95,
  type = "perc",
  R = 1000,
  histogram = FALSE,
  digits = 3,
  reportIncomplete = FALSE
)
}
\arguments{
\item{x}{Either a two-way table or a two-way matrix.
Can also be a vector of observations of an ordinal variable.}

\item{g}{If \code{x} is a vector, \code{g} is the vector of observations for
the grouping, nominal variable.}

\item{group}{If \code{x} is a table or matrix, \code{group} indicates whether
the \code{"row"} or the \code{"column"} variable is
the nominal, grouping variable.}

\item{verbose}{If \code{TRUE}, prints statistics for each
comparison.}

\item{progress}{If \code{TRUE}, prints a message as each comparison is
conducted.}

\item{ci}{If \code{TRUE}, returns confidence intervals by bootstrap.
May be slow.}

\item{conf}{The level for the confidence interval.}

\item{type}{The type of confidence interval to use.
Can be any of "\code{norm}", "\code{basic}", 
              "\code{perc}", or "\code{bca}".
Passed to \code{boot.ci}.}

\item{R}{The number of replications to use for bootstrap.}

\item{histogram}{If \code{TRUE}, produces a histogram of bootstrapped values.}

\item{digits}{The number of significant digits in the output.}

\item{reportIncomplete}{If \code{FALSE} (the default),
\code{NA} will be reported in cases where there
are instances of the calculation of the statistic
failing during the bootstrap procedure.}
}
\value{
A single statistic, Freeman's theta.
        Or a small data frame consisting of Freeman's theta,
        and the lower and upper confidence limits.
}
\description{
Calculates Freeman's theta for a table with one ordinal
             variable and one nominal variable; confidence intervals
             by bootstrap.
}
\details{
Freeman's coefficent of differentiation (theta)
          is used as a measure of association
          for a two-way
          table with one ordinal and one nominal variable.
          See Freeman (1965).
          
          Currently, the function makes no provisions for \code{NA}
          values in the data.  It is recommended that \code{NA}s be removed
          beforehand.
          
          Because theta is always positive, if \code{type="perc"},
          the confidence interval will
          never cross zero, and should not
          be used for statistical inference.
          However, if \code{type="norm"}, the confidence interval
          may cross zero.

          When theta is close to 0 or very large,
          or with small counts in some cells,
          the confidence intervals 
          determined by this
          method may not be reliable, or the procedure may fail.
}
\examples{
data(Breakfast)
library(coin)
chisq_test(Breakfast, scores = list("Breakfast" = c(-2, -1, 0, 1, 2)))
freemanTheta(Breakfast)

### Example from Freeman (1965), Table 10.6
Counts = c(1, 2, 5, 2, 0, 10, 5, 5, 0, 0, 0, 0, 2, 2, 1, 0, 0, 0, 2, 3)
Matrix = matrix(Counts, byrow=TRUE, ncol=5,
                dimnames = list(Marital.status = c("Single", "Married",
                                                   "Widowed", "Divorced"),
                                Social.adjustment = c("5","4","3","2","1")))
Matrix
freemanTheta(Matrix)

### Example after Kruskal Wallis test
data(PoohPiglet)
kruskal.test(Likert ~ Speaker, data = PoohPiglet)
freemanTheta(x = PoohPiglet$Likert, g = PoohPiglet$Speaker)

### Same data, as table of counts
data(PoohPiglet)
XT = xtabs( ~ Speaker + Likert , data = PoohPiglet)
freemanTheta(XT)

### Example from Freeman (1965), Table 10.7
Counts = c(52, 28, 40, 34, 7, 9, 16, 10, 8, 4, 10, 9, 12,6, 7, 5)
Matrix = matrix(Counts, byrow=TRUE, ncol=4,
                dimnames = list(Preferred.trait = c("Companionability",
                                                    "PhysicalAppearance",
                                                    "SocialGrace",
                                                    "Intelligence"),
                                Family.income = c("4", "3", "2", "1")))
Matrix
freemanTheta(Matrix, verbose=TRUE)

}
\references{
Freeman, L.C. 1965. Elementary Applied Statistics for Students
            in Behavioral Science. Wiley.

            \url{https://rcompanion.org/handbook/H_11.html}
}
\seealso{
\code{\link{epsilonSquared}}
}
\author{
Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
}
\concept{Cochran-Armitage}
\concept{Freeman's theta}
\concept{Kruskal-Wallis}
