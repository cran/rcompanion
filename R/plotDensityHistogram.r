#' @title Histogram with a density curve
#'
#' @description Produces a histogram for a vector of values and adds a density
#'              curve of the distribution. 
#' 
#' @param x A vector of values.
#' @param prob If \code{FALSE}, then counts are displayed in the histogram.
#'        If \code{TRUE}, then the density is shown.
#' @param col The color of the histogram bars.
#' @param main The title displayed for the plot.
#' @param linecol The color of the line in the plot.
#' @param lwd The width of the line in the plot.
#' @param adjust Passed to \code{\link{density}}. A lower value makes the density
#'               plot smoother.
#' @param bw Passed to \code{\link{density}}.
#' @param kernel Passed to \code{\link{density}}.
#' @param ... Other arguments passed to \code{\link{hist}}.
#' 
#' @details  The function relies on the \code{hist} function. The density curve
#'           relies on the \code{density} function.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/C_04.html}
#' @seealso \code{\link{plotNormalHistogram}} \code{\link{plotNormalDensity}}
#' @concept density distribution histogram
#' @return Produces a plot. Returns nothing.
#'         
#' @examples
#' ### Plot of residuals from a model fit with lm
#' data(Catbus)
#' model = lm(Steps ~ Sex + Teacher,
#'            data = Catbus)
#' plotDensityHistogram(residuals(model))          
#' 
#' @importFrom stats lm density
#' @importFrom graphics hist lines
#' 
#' @export

plotDensityHistogram = 
   function(x, prob=FALSE, col="gray",
                           main="", linecol="black", lwd=2, 
                           adjust=1, bw="nrd0",
                           kernel="gaussian", ...)
   {
   x = x[!is.na(x)]
   Hissy = hist(x, plot=FALSE, ...)
   if(prob==TRUE){
      Range = density(x, adjust=adjust, bw=bw, kernel=kernel)$x
      Dens  = density(x, adjust=adjust, bw=bw, kernel=kernel)$y
      Max   = max(c(Hissy$density, max(Dens)))
   }
   if(prob==FALSE){
      Factor = mean(Hissy$counts / Hissy$density)
      Range = density(x, adjust=adjust, bw=bw, kernel=kernel)$x
      Dens = density(x, adjust=adjust, bw=bw, kernel=kernel)$y * Factor
      Max = max(c(Hissy$count, max(Dens)))
   }
   hist(x,
        prob = prob,
        col  = col, 
        main = main,
        ylim = c(0, Max),
        ...)
   lines(Range, Dens,
         col = linecol,
         lwd = lwd)
}
