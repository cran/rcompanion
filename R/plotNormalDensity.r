#' @title Density plot with a normal curve
#'
#' @description Produces a density plot for a vector of values and adds 
#'              a normal curve with the same mean and standard deviation. The 
#'              plot can be used to quickly compare the distribution of data
#'              to a normal distribution. 
#' 
#' @param x A vector of values.
#' @param col1 The color of the density plot. Usually not visible.
#' @param col2 The color of the density polygon.
#' @param col3 The color of the normal line.
#' @param border The color of the border around the density polygon.
#' @param main The title displayed for the plot.
#' @param lwd The width of the line in the plot.
#' @param length The number of points in the line in the plot.
#' @param adjust Passed to  \code{\link{density}}. 
#'        A lower value makes the density plot smoother.
#' @param bw Passed to  \code{\link{density}}.
#' @param kernel Passed to  \code{\link{density}}.
#' @param ... Other arguments passed to \code{\link{plot}}.
#' 
#' @details  The function plots a polygon based on the \code{density} function. 
#'           The normal curve has the same mean and standard deviation as the 
#'           values in the vector.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/I_01.html}
#' @seealso \code{\link{plotNormalHistogram}} \code{\link{plotDensityHistogram}}
#' @concept normal distribution density
#' @return Produces a plot. Returns nothing.
#'         
#' @examples
#' ### Plot of residuals from a model fit with lm
#' data(Catbus)
#' model = lm(Steps ~ Gender + Teacher,
#'            data = Catbus)
#'  plotNormalDensity(residuals(model))          
#' 
#' @importFrom stats sd dnorm density
#' @importFrom graphics plot lines polygon
#'
#' @export

plotNormalDensity = 
   function(x, col1="white", col2="gray", col3="blue",
                           border=NA, main="", lwd=2, length=1000,
                           adjust=1, bw="nrd0",
                           kernel="gaussian", ...)
   {
   x = x[!is.na(x)]
   Range = seq(min(x), max(x), length = length)
   Norm  = dnorm(Range, mean = mean(x), sd   = sd(x))
   Dens = density(x, adjust=adjust, bw=bw, kernel=kernel)
   Max = max(c(Norm, Dens$y))
   plot(Dens, col=col1, ylim=c(0, Max), main=main, ...)
   polygon(Dens, col=col2, border=border)
   lines(Range, Norm, col = col3, lwd = lwd)
}
