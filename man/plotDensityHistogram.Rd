% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotDensityHistogram.r
\name{plotDensityHistogram}
\alias{plotDensityHistogram}
\title{Histogram with a density curve}
\usage{
plotDensityHistogram(
  x,
  prob = FALSE,
  col = "gray",
  main = "",
  linecol = "black",
  lwd = 2,
  adjust = 1,
  bw = "nrd0",
  kernel = "gaussian",
  ...
)
}
\arguments{
\item{x}{A vector of values.}

\item{prob}{If \code{FALSE}, then counts are displayed in the histogram.
If \code{TRUE}, then the density is shown.}

\item{col}{The color of the histogram bars.}

\item{main}{The title displayed for the plot.}

\item{linecol}{The color of the line in the plot.}

\item{lwd}{The width of the line in the plot.}

\item{adjust}{Passed to \code{\link{density}}. A lower value makes the density
plot smoother.}

\item{bw}{Passed to \code{\link{density}}.}

\item{kernel}{Passed to \code{\link{density}}.}

\item{...}{Other arguments passed to \code{\link{hist}}.}
}
\value{
Produces a plot. Returns nothing.
}
\description{
Produces a histogram for a vector of values and adds a density
             curve of the distribution.
}
\details{
The function relies on the \code{hist} function. The density curve
          relies on the \code{density} function.
}
\examples{
### Plot of residuals from a model fit with lm
data(Catbus)
model = lm(Steps ~ Gender + Teacher,
           data = Catbus)
plotDensityHistogram(residuals(model))          

}
\references{
\url{https://rcompanion.org/handbook/C_04.html}
}
\seealso{
\code{\link{plotNormalHistogram}}, 
         \code{\link{plotNormalDensity}}
}
\author{
Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
}
\concept{density plot}
\concept{histogram}
\concept{normal distribution}
