#' @title Histogram with a normal curve
#'
#' @description Produces a histogram for a vector of values and adds a normal
#'              curve with the same mean and standard deviation. The plot
#'              can be used to quickly compare the distribution of data
#'              to a normal distribution. 
#' 
#' @param x A vector of values.
#' @param prob If \code{FALSE}, then counts are displayed in the histogram.
#'        If \code{TRUE}, then the density is shown.
#' @param col The color of the histogram bars.
#' @param main The title displayed for the plot.
#' @param linecol The color of the line in the plot.
#' @param lwd The width of the line in the plot.
#' @param length The number of points in the line in the plot.
#' @param ... Other arguments passed to \code{\link{hist}}.
#' 
#' @details  The function relies on the \code{hist} function. The normal curve
#'           has the same mean and standard deviation as the values in the
#'           vector.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/I_01.html}
#' @seealso \code{\link{plotNormalDensity}} \code{\link{plotDensityHistogram}}
#' @concept normal distribution histogram
#' @return Produces a plot. Returns nothing.
#'         
#' @examples
#' ### Plot of residuals from a model fit with lm
#' data(Catbus)
#' model = lm(Steps ~ Sex + Teacher,
#'            data = Catbus)
#'  plotNormalHistogram(residuals(model))          
#' 
#' @importFrom stats sd dnorm
#' @importFrom graphics hist lines
#' 
#' @export

plotNormalHistogram = 
   function(x, prob=FALSE, col="gray",
                            main="", linecol="blue", lwd=2, 
                            length = 1000, ...)
   {
   x = x[!is.na(x)]
   Range = seq(min(x), max(x), length = length)
   Hissy = suppressWarnings(hist(x, plot=FALSE, ...))
   if(prob==TRUE){
      Norm  = dnorm(Range, mean = mean(x), sd   = sd(x))
      Max = max(c(Hissy$density, max(Norm)))
   }
   if(prob==FALSE){
      Factor = mean(Hissy$counts[1] / Hissy$density[1])
      Norm  = dnorm(Range, mean = mean(x), sd   = sd(x)) * Factor
      Max = max(c(Hissy$count, max(Norm)))
   }
   hist(x,
        prob = prob,
        col  = col, 
        main = main,
        ylim = c(0, Max),
        ...)
   lines(Range, Norm,
         col = linecol,
         lwd = lwd)
}