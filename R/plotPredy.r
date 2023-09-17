#' @title Plot a predicted line from a bivariate model
#'
#' @description Plots the best fit line for a model with one y variable and
#'              one x variable, or with one y variable and polynomial
#'              x variables.
#' 
#' @param data   The name of the data frame.
#' @param x      The name of the x variable.
#' @param y      The name of the y variable.
#' @param model  The name of the model object.
#' @param order  If plotting a polynomial function, the order of the polynomial.
#'               Otherwise can be left as \code{1}.
#' @param x2     If applicable, the name of the second order 
#'               polynomial x variable.
#' @param x3     If applicable, the name of the third order 
#'               polynomial x variable.
#' @param x4     If applicable, the name of the fourth order 
#'               polynomial x variable.
#' @param x5     If applicable, the name of the fifth order 
#'               polynomial x variable.
#' @param pch    The shape of the plotted data points.
#' @param xlab   The label for the x-axis.
#' @param ylab   The label for the y-axis.
#' @param length The number of points used to draw the line.
#' @param lty    The style of the plotted line.
#' @param lwd    The width of the plotted line.
#' @param col    The col of the plotted line.
#' @param type   Passed to \code{predict}. Required for certain models.
#' @param ...    Other arguments passed to \code{plot}.
#' 
#' @details  Any model for which \code{predict()} is defined can be used.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{https://rcompanion.org/handbook/I_10.html}
#' 
#' @return Produces a plot. Returns nothing.
#'         
#' @examples
#' ### Plot of linear model fit with lm
#' data(BrendonSmall)
#' model = lm(Weight ~ Calories, data = BrendonSmall) 
#' plotPredy(data  = BrendonSmall,
#'           y     = Weight,
#'           x     = Calories,
#'           model = model,
#'           xlab  = "Calories per day",
#'           ylab  = "Weight in kilograms")
#'            
#' ### Plot of polynomial model fit with lm
#' data(BrendonSmall)
#' BrendonSmall$Calories2 = BrendonSmall$Calories * BrendonSmall$Calories
#' model = lm(Sodium ~ Calories + Calories2, data = BrendonSmall) 
#' plotPredy(data  = BrendonSmall,
#'           y     = Sodium,
#'           x     = Calories,
#'           x2    = Calories2,
#'           model = model,
#'           order = 2,
#'           xlab  = "Calories per day",
#'           ylab  = "Sodium intake per day")
#' 
#' ### Plot of quadratic plateau model fit with nls
#' data(BrendonSmall)
#' quadplat = function(x, a, b, clx) {
#'           ifelse(x  < clx, a + b * x   + (-0.5*b/clx) * x   * x,
#'                            a + b * clx + (-0.5*b/clx) * clx * clx)}
#' model = nls(Sodium ~ quadplat(Calories, a, b, clx),
#'             data = BrendonSmall,
#'             start = list(a   = 519,
#'                          b   = 0.359,
#'                          clx = 2304))
#' plotPredy(data  = BrendonSmall,
#'           y     = Sodium,
#'           x     = Calories,
#'           model = model,
#'           xlab  = "Calories per day",
#'           ylab  = "Sodium intake per day")
#'
#' ### Logistic regression example requires type option
#' data(BullyHill)
#' Trials = cbind(BullyHill$Pass, BullyHill$Fail)
#' model.log = glm(Trials ~ Grade, data = BullyHill,
#'                 family = binomial(link="logit"))
#' plotPredy(data  = BullyHill,
#'           y     = Percent,
#'           x     = Grade,
#'           model = model.log,
#'           type  = "response",
#'           xlab  = "Grade",
#'           ylab  = "Proportion passing")
#' 
#' @importFrom graphics plot lines
#' @importFrom stats lm predict nls formula
#' 
#' @export

plotPredy = function (data, x, y, model, order=1,
                       x2=NULL, x3=NULL, x4=NULL, x5=NULL,
                       pch=16, xlab="X", ylab="Y", length = 1000,
                       lty=1, lwd=2, col="blue", type=NULL, ...) {
   x=as.character(substitute(x))
   y=as.character(substitute(y))
   if(order>1){x2=as.character(substitute(x2))}
   if(order>2){x3=as.character(substitute(x3))}
   if(order>3){x4=as.character(substitute(x4))}
   if(order>4){x5=as.character(substitute(x5))}
   Formula = formula(paste(y, "~", x, sep=" "))
   i = as.numeric(seq(min(data[,x]), max(data[,x]), len=length))
   if(order==1){D1 = data.frame(X=i)
                names(D1)[1] = x}
   if(order==2){D1 = data.frame(X=i, X2=i*i)
                names(D1)[1] = x
                names(D1)[2] = x2}
   if(order==3){D1 = data.frame(X=i, X2=i*i, X3=i*i*i)
                names(D1)[1] = x
                names(D1)[2] = x2
                names(D1)[3] = x3}
   if(order==4){D1 = data.frame(X=i, X2=i*i, X3=i*i*i, X4=i*i*i*i)
   names(D1)[1] = x
   names(D1)[2] = x2
   names(D1)[3] = x3
   names(D1)[4] = x4}
   if(order==5){D1 = data.frame(X=i, X2=i*i, X3=i*i*i, X4=i*i*i*i, X5=i*i*i*i*i)
   names(D1)[1] = x
   names(D1)[2] = x2
   names(D1)[3] = x3
   names(D1)[4] = x4
   names(D1)[5] = x5}
   if(is.null(type)) {predy = predict(model, D1)}
   if(!is.null(type)) {predy = predict(model, D1, type=type)}
   plot(Formula, data = data, pch=pch, xlab=xlab, ylab=ylab, ...)
   lines(i, predy,
         lty=lty,                                  
         lwd=lwd, 
         col=col)
}