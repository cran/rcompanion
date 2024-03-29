% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wilcoxonZ.r
\name{wilcoxonZ}
\alias{wilcoxonZ}
\title{Wilcoxon z statistic}
\usage{
wilcoxonZ(
  x,
  y = NULL,
  mu = 0,
  paired = FALSE,
  exact = FALSE,
  correct = FALSE,
  digits = 3
)
}
\arguments{
\item{x}{A vector of observations.}

\item{y}{For the two-sample and paired cases,
a second vector of observations.}

\item{mu}{For the one-sample case,
the value to compare \code{x} to, as in \code{wilcox.test}}

\item{paired}{As used in \code{wilcox.test}.}

\item{exact}{As used in \code{wilcox.test}, 
default here is \code{FALSE}.}

\item{correct}{As used in \code{wilcox.test}, 
default here is \code{FALSE}.}

\item{digits}{The number of significant digits in the output.}
}
\value{
A single statistic, \code{z}.
}
\description{
Calculates the z statistic for a Wilcoxon
             two-sample, paired, or one-sample test.
}
\details{
This function uses code from \code{wilcox.test},
          and reports the \code{z} statistic,
          which is calculated by the original function
          but isn't returned.
          
          The returned value will be NA if the function attempts an
          exact test.
          
          For the paired case, the observations in \code{x} and
          and \code{y} should be ordered such that the
          first observation in \code{x} is paired with the first observation
          in \code{y}, and so on.
}
\examples{
data(Pooh)
wilcoxonZ(x = Pooh$Likert[Pooh$Time==1], y = Pooh$Likert[Pooh$Time==2],
          paired=TRUE, exact=FALSE, correct=FALSE)

}
\author{
Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}, 
        R Core Team
}
\concept{Wilcoxon signed rank}
\concept{Wilcoxon-Mann-Whitney}
