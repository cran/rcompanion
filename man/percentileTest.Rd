% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/percentileTest.r
\name{percentileTest}
\alias{percentileTest}
\title{Test of percentiles by permutation test}
\usage{
percentileTest(
  formula = NULL,
  data = NULL,
  x = NULL,
  y = NULL,
  test = "median",
  tau = 0.5,
  type = 7,
  threshold = NA,
  comparison = "<",
  r = 1000,
  digits = 4,
  progress = "TRUE"
)
}
\arguments{
\item{formula}{A formula indicating the response variable and
the independent variable. e.g. y ~ group.}

\item{data}{The data frame to use.}

\item{x}{If no formula is given, the response variable for one group.}

\item{y}{The response variable for the other group.}

\item{test}{The statistic to compare between groups.  Can be
\code{"median"}, \code{"percentile"}, \code{"iqr"},
\code{"proportion"},
\code{"mean"}, or \code{"variance"}.}

\item{tau}{If \code{"percentile"} is chosen as the \code{test},
\code{tau} indicates the percentile to test.  Expressed
as a quantile.  That is, 0.5 indicates a test for medians.
0.75 indicates a test for 75th percentiles.}

\item{type}{The \code{type} value passed to the \code{quantile} function.}

\item{threshold}{If \code{"proportion"} is chosen as the \code{test},
\code{threshold} indicates the value of the dependent variable
to use as the threshold.  For example, to test if there is a 
different in the proportion of observations below $10,000,
\code{threshold = 10000} would be used.}

\item{comparison}{If \code{"proportion"} is chosen as the \code{test},
\code{comparison} indicates the inequality to use.  Options are
\code{"<"}, \code{"<="}, \code{">"}, \code{">="}, or , \code{"=="}}

\item{r}{The number of replicates in the permutation test.}

\item{digits}{The number of significant digits in the output.}

\item{progress}{If \code{TRUE}, prints a dot for every 1 percent of 
progress while conducting the test.}
}
\value{
A list of three data frames with
        the data used, a summary for each group,
        and the p-value from the test.
}
\description{
Conducts a permutation test to compare two groups for medians,
             percentiles, or proportion below a threshold value.
}
\details{
The function will test for a difference in medians, percentiles,
         interquartile ranges, proportion of observations above or below
         some threshold value, 
         means, or variances between two groups
         by permutation test.
         
         The permutation test simply permutes the observed values over the
         two groups and counts how often the calculated statistic is
         at least as extreme as the original observed statistic.
         
         The input should include either \code{formula} and \code{data};
         or \code{x} and \code{y}.
         
         The function removes cases with NA in any of the variables.
         
         If the independent variable has more than two groups,
         only the first two levels of the factor variable will be used.
         
         The p-value returned is a two-sided test.
}
\note{
The parsing of the formula is simplistic. 
         The first variable on the
         left side is used as the measurement variable.  
         The first variable on the
         right side is used for the independent variable.
}
\examples{
data(BrendonSmall)
percentileTest(Sodium ~ Instructor, 
               data=BrendonSmall, 
               test="median")

percentileTest(Sodium ~ Instructor, 
               data=BrendonSmall, 
               test="percentile", 
               tau = 0.75)

percentileTest(Sodium ~ Instructor, 
               data=BrendonSmall, 
               test="proportion", 
               threshold = 1300)
               
}
\references{
\url{https://rcompanion.org/handbook/F_15.html}
}
\author{
Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
}
\concept{median}
\concept{percentile}
\concept{permutation test}
\concept{quantile}
