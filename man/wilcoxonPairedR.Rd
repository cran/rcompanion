% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wilcoxonPairedR.r
\name{wilcoxonPairedR}
\alias{wilcoxonPairedR}
\title{r effect size for Wilcoxon two-sample paired signed-rank test}
\usage{
wilcoxonPairedR(
  x,
  g = NULL,
  adjustn = TRUE,
  coin = FALSE,
  ci = FALSE,
  conf = 0.95,
  type = "perc",
  R = 1000,
  histogram = FALSE,
  cases = TRUE,
  digits = 3,
  ...
)
}
\arguments{
\item{x}{A vector of observations.}

\item{g}{The vector of observations for
the grouping, nominal variable.
Only the first two levels of the nominal variable are used.
The data must be ordered so that the first observation of the
of the first group is paired with the first observation of the
second group.}

\item{adjustn}{If \code{TRUE}, reduces the sample size in the calculation
of \code{r} by the number of tied pairs.}

\item{coin}{If \code{FALSE}, the default, the Z value
is extracted from a function similar to the
\code{wilcox.test} function in the stats package.
If \code{TRUE}, the Z value
is extracted from the \code{wilcox_test} function in the
coin package.  This method may be much slower, especially
if a confidence interval is produced.}

\item{ci}{If \code{TRUE}, returns confidence intervals by bootstrap.
May be slow.}

\item{conf}{The level for the confidence interval.}

\item{type}{The type of confidence interval to use.
Can be any of "\code{norm}", "\code{basic}", 
              "\code{perc}", or "\code{bca}".
Passed to \code{boot.ci}.}

\item{R}{The number of replications to use for bootstrap.}

\item{histogram}{If \code{TRUE}, produces a histogram of bootstrapped values.}

\item{cases}{By default the \code{N} used in the formula for \code{r}
is the number of pairs.  If \code{cases=FALSE},
the \code{N} used in the formula for \code{r}
is the total number of observations, as some sources suggest.}

\item{digits}{The number of significant digits in the output.}

\item{...}{Additional arguments passed to the \code{wilcoxsign_test} 
function.}
}
\value{
A single statistic, r.
        Or a small data frame consisting of r,
        and the lower and upper confidence limits.
}
\description{
Calculates r effect size
             for a Wilcoxon two-sample paired signed-rank test; 
             confidence intervals by bootstrap.
}
\details{
r  is calculated as Z divided by 
          square root of the number of observations in one group. This 
          results in a statistic that ranges from -1 to 1.
          This range doesn't hold if \code{cases=FALSE}.
          
          This statistic typically reports a smaller effect size
          (in absolute value) than does
          the matched-pairs rank biserial correlation coefficient 
          (\code{wilcoxonPairedRC}), and may not reach a value
          of -1 or 1 if there are ties in the paired differences.

          Currently, the function makes no provisions for \code{NA}
          values in the data.  It is recommended that \code{NA}s be removed
          beforehand.

          When the data in the first group are greater than
          in the second group, r is positive.
          When the data in the second group are greater than
          in the first group, r is negative.
          Be cautious with this interpretation, as R will alphabetize
          groups if \code{g} is not already a factor.
          
          When r is close to extremes,
          or with small counts in some cells,
          the confidence intervals 
          determined by this
          method may not be reliable, or the procedure may fail.
}
\section{Acknowledgments}{

         My thanks to
         Peter Stikker for the suggestion to adjust the sample size for ties.
}

\examples{
data(Pooh)
Time1 = Pooh$Likert[Pooh$Time==1]
Time2 = Pooh$Likert[Pooh$Time==2]
wilcox.test(x = Time1, y = Time2, paired=TRUE, exact=FALSE)
wilcoxonPairedR(x = Pooh$Likert, g = Pooh$Time)

}
\references{
\url{https://rcompanion.org/handbook/F_06.html}
}
\seealso{
\code{\link{wilcoxonPairedRC}}
}
\author{
Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
}
\concept{Wilcoxon signed rank}
\concept{confidence interval}
\concept{effect size}
