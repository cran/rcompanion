#' @title Cate-Nelson models for bivariate data with a fixed critical Y value
#'
#' @description Produces critical-x values for bivariate data
#'              according to a Cate-Nelson analysis
#'              for a given critical Y value.
#' 
#' @param x A vector of values for the x variable.
#' @param y A vector of values for the y variable.
#' @param cly = Critical Y value.
#' @param plotit If \code{TRUE}, produces plots of the output.
#' @param hollow If \code{TRUE}, uses hollow circles on the plot to indicate
#'               data not fitting the model.
#' @param xlab The label for the x-axis.
#' @param ylab The label for the y-axis.
#' @param trend \code{"postive"} if the trend of y vs. x is generally
#'              positive. \code{"negative"} if negative.
#' @param clx Indicates which of the listed critical x values 
#'            should be chosen as the critical x value for the plot.
#' @param outlength Indicates the number of potential critical x values
#'                   to display in the output.
#' @param sortstat The statistic to sort by.  Any of \code{"error"} 
#'                 (the default),
#'                 \code{"phi"}, \code{"fisher"}, or \code{"pearson"}.                       
#'             
#' @details  Cate-Nelson analysis divides bivariate data into two groups.
#'           For data with a positive trend, one group has a
#'           large \code{x} value associated with a large \code{y} value, and   
#'           the other group has a small \code{x} value associated with a small   
#'           \code{y} value. For a negative trend, a small \code{x} is
#'           associated with a large \code{y}, and so on. 
#'           
#'           The analysis is useful for bivariate data which don't conform well
#'           to linear, curvilinear, or plateau models.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{http://rcompanion.org/rcompanion/h_02.html}
#' 
#' @seealso \code{\link{cateNelson}}
#' 
#' @concept Cate-Nelson
#' @concept agronomy
#' 
#' @return A data frame of statistics from the analysis:
#'         critical level for x, critical value for y, the
#'         number of observations in each of the quadrants (I, II, III, IV),
#'         the number of observations that conform with the model, 
#'         the number of observations that do not conform to the model,
#'         the proportion of observations that conform with the model,
#'         the proportion of observations that do not conform to the model,
#'         a p-value for the Fisher exact test for the data divided into
#'         the groups indicated by the model,
#'         phi for the data divided into
#'         the groups indicated by the model,
#'         and Pearson's chi-square for the data divided into
#'         the groups indicated by the model.
#'                   
#' @examples
#' data(Nurseries)
#' cateNelsonFixedY(x          = Nurseries$Size,
#'                  y          = Nurseries$Proportion,
#'                  cly        = 0.70,
#'                  plotit     = TRUE,
#'                  hollow     = TRUE,
#'                  xlab       = "Nursery size in hectares",
#'                  ylab       = "Proportion of good practices adopted",
#'                  trend      = "positive",
#'                  clx        = 1,
#'                  outlength  = 15)
#' 
#' @importFrom graphics plot abline text
#' @importFrom stats fisher.test chisq.test
#' @importFrom utils head
#' 
#' @export
#' 
cateNelsonFixedY = 
  function(x, y, cly=0.95, plotit=TRUE, hollow=TRUE, xlab="X", ylab="Y", 
           trend="positive", clx=1,
           outlength=20, sortstat="error")
  {
    N = length(x)
    n = N-1
    XS = sort(x)
    
    Out = data.frame(Critx  = as.numeric(rep(0.00, n)),
                     Crity  = as.numeric(rep(0.00, n)),
                     Q1     = as.integer(rep(0, n)),
                     Q2     = as.integer(rep(0, n)),
                     Q3     = as.integer(rep(0, n)),
                     Q4     = as.integer(rep(0, n)),
                     Model  = as.integer(rep(0, n)),
                     Error  = as.integer(rep(0, n)),
                     N      = as.integer(rep(0, n)),
                     pQ1    = as.numeric(rep(0, n)),
                     pQ2    = as.numeric(rep(0, n)),
                     pQ3    = as.numeric(rep(0, n)),
                     pQ4    = as.numeric(rep(0, n)),
                     pModel = as.numeric(rep(0, n)),
                     pError = as.numeric(rep(0, n)),
                     Fisher.p = as.numeric(rep(0, n)),
                     Pearson.chisq = as.numeric(rep(0, n)),
                     Pearson.p = as.numeric(rep(0, n)),
                     phi = as.numeric(rep(0, n)))
    
    for(i in c(1:(N-1))){
      Out[i, 1] = mean(XS[i+1], XS[i])
      Out[i, 2] = cly
      Out[i, 3] = sum(y>Out$Crity[i]  & x<Out$Critx[i])   # Q1
      Out[i, 4] = sum(y>=Out$Crity[i] & x>=Out$Critx[i])  # Q2
      Out[i, 5] = sum(y<Out$Crity[i]  & x>Out$Critx[i])   # Q3
      Out[i, 6] = sum(y<=Out$Crity[i] & x<=Out$Critx[i])  # Q4
      if (trend=="positive")
      {
      Out[i, 7] = sum(c(Out$Q2[i], Out$Q4[i]))            # Model
      Out[i, 8] = sum(c(Out$Q1[i], Out$Q3[i]))            # Error
      }
      if (trend=="negative")
      {
        Out[i, 7] = sum(c(Out$Q1[i], Out$Q3[i]))            # Model
        Out[i, 8] = sum(c(Out$Q2[i], Out$Q4[i]))            # Error
      }
      Out[i, 9] = sum(c(Out$Q1[i], Out$Q2[i], 
                        Out$Q3[i], Out$Q4[i]))            # N
      Out[i, 10]  = round(Out$Q1[i] / Out$N[i], 3)        # percents
      Out[i, 11]  = round(Out$Q2[i] / Out$N[i], 3)        # percents
      Out[i, 12]  = round(Out$Q3[i] / Out$N[i], 3)        # percents
      Out[i, 13]  = round(Out$Q4[i] / Out$N[i], 3)        # percents  
      Out[i, 14]  = round(Out$Model[i] / Out$N[i], 3)     # percents
      Out[i, 15]  = round(Out$Error[i] / Out$N[i], 3)     # percents 
      
      M=matrix(c(Out$Q1[i],Out$Q2[i],Out$Q4[i],Out$Q3[i]), nrow=2, byrow=TRUE)
      Out[i, 16]  = signif(fisher.test(M)$p.value, 4)
      Out[i, 17]  = signif(suppressWarnings(chisq.test(M,correct=TRUE)$statistic), 4)
      Out[i, 18]  = signif(suppressWarnings(chisq.test(M,correct=TRUE)$p.value), 4)
      
      MM = M / sum(M)
      a=MM[1,1]; b=MM[1,2]; c=MM[2,1]; d=MM[2,2]
      Out[i, 19] = round((a- (a+b)*(a+c))/sqrt((a+b)*(c+d)*(a+c)*(b+d)),3)
    }

    Out2 = Out[order(Out$Error),]
    if(sortstat=="phi"){Out2 = Out[order(abs(Out$phi), decreasing=TRUE),]}
    if(sortstat=="model"){Out2 = Out[order(Out$Model, decreasing=TRUE),]}
    if(sortstat=="fisher"){Out2 = Out[order(Out$Fisher.p),]}
    if(sortstat=="pearson"){Out2 = Out[order(Out$Pearson.chisq, decreasing=TRUE),]}
    row.names(Out2) <- 1:nrow(Out2)
    
  if(plotit){
    
    plot(pModel ~ Critx,
         data=Out,
         xlab="Critical-x value",
         ylab="Proportion of obs fitting model")
    
    plot(abs(phi) ~ Critx,
         data=Out,
         xlab="Critical-x value",
         ylab="abs(phi)")
    
    plot(Pearson.chisq ~ Critx,
         data=Out,
         xlab="Critical-x value",
         ylab="Pearson chi-square")
    
    CLX = Out2$Critx[clx]
    CLY = Out2$Crity[clx]
    
    if (trend=="positive")
    {
      pchi=as.integer(rep(16,n))
      pchi[(x<CLX)&(y<CLY)] = 16
      pchi[(x>CLX)&(y>CLY)] = 16
      pchi[(x>CLX)&(y<CLY)] = 1
      pchi[(x<CLX)&(y>CLY)] = 1
    }
    if (trend=="negative")
    {
      pchi=as.integer(rep(16,n))
      pchi[(x<CLX)&(y<CLY)] = 1
      pchi[(x>CLX)&(y>CLY)] = 1
      pchi[(x>CLX)&(y<CLY)] = 16
      pchi[(x<CLX)&(y>CLY)] = 16
    }
    if (hollow==FALSE)
    {
      pchi=as.integer(rep(16,n))
    }
    plot(x, y, pch=pchi,
         xlab=xlab,
         ylab=ylab)
    abline(v=CLX, col="blue")
    abline(h=CLY, col="blue")
  }
    
  return(head(Out2, outlength))
}