#' @title Cate-Nelson models for bivariate data
#'
#' @description Produces critical-x and critical-y values for bivariate data
#'              according to a Cate-Nelson analysis.
#' 
#' @param x A vector of values for the x variable.
#' @param y A vector of values for the y variable.
#' @param plotit If \code{TRUE}, produces plots of the output.
#' @param hollow If \code{TRUE}, uses hollow circles on the plot to indicate
#'               data not fitting the model.
#' @param xlab The label for the x-axis.
#' @param ylab The label for the y-axis.
#' @param trend \code{"postive"} if the trend of y vs. x is generally
#'              positive. \code{"negative"} if negative.
#' @param clx Indicates which of the listed critical x values 
#'            should be chosen as the critical x value for the final model.
#' @param cly Indicates which of the listed critical y values 
#'            should be chosen as the critical y value for the final model.
#' @param xthreshold Indicates the proportion of potential critical x values
#'                   to display in the output. A value of \code{1} would display
#'                   all of them.
#' @param ythreshold Indicates the proportion of potential critical y values
#'                   to display in the output. A value of \code{1} would display
#'                   all of them.
#' @param progress If \code{TRUE}, prints an indicator of progress as for loops
#'                 progress.
#' @param verbose If \code{FALSE}, suppresses printed output of tables.
#' @param listout If \code{TRUE}, outputs a list of data frames instead of a
#'                a single data frame.  This allows a data frame of critical
#'                values and associated statistics to be extracted, for example
#'                if one would want to sort by Cramer's V.                             
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
#'           This function will fail if either of the largest two or smallest 
#'           two x values are identical.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/rcompanion/h_02.html}
#' @concept Cate Nelson bivariate soil agronomy
#' @return A data frame of statistics from the analysis: number of observations,
#'         critical level for x, sum of squares, critical value for y, the
#'         number of observations in each of the quadrants (I, II, III, IV),
#'         the number of observations that conform with the model, 
#'         the proportion of observations that conform with the model, 
#'         the number of observations that do not conform to the model,
#'         the proportion of observations that do not conform to the model,
#'         a p-value for the Fisher exact test for the data divided into
#'         the groups indicated by the model,
#'         and Cramer's V for the data divided into
#'         the groups indicated by the model.
#'         
#'         Output also includes printed lists of critical values,
#'         explanation of the values in the data frame,
#'         and plots: y vs. x; sum of squares vs. critical x value;
#'         the number of observations that do not conform to the model vs.
#'         critical y value; 
#'         and y vs. x with the critical values shown as lines on the plot, 
#'         and the quadrants labeled.
#' 
#' @note  The method in this function follows \cite{Cate, R. B., & Nelson, L.A. 
#'        (1971). A simple statistical procedure for partitioning soil test 
#'        correlation data into two classes. Soil Science Society of America 
#'        Proceedings 35, 658-660.}
#'           
#'        An earlier version of this function was published in \cite{Mangiafico, 
#'        S.S. 2013. Cate-Nelson Analysis for Bivariate Data Using 
#'        R-project. J.of Extension 51:5, 5TOT1.}
#'                   
#' @examples
#' data(Nurseries)
#' cateNelson(x          = Nurseries$Size,
#'            y          = Nurseries$Proportion,
#'            plotit     = TRUE,
#'            hollow     = TRUE,
#'            xlab       = "Nursery size in hectares",
#'            ylab       = "Proportion of good practices adopted",
#'            trend      = "positive",
#'            clx        = 1,
#             cly        = 1,
#'            xthreshold = 0.10,
#'            ythreshold = 0.15)
#' 
#' @importFrom graphics plot abline text
#' @importFrom stats anova fisher.test lm
#' 
#' @export
#' 
cateNelson = 
  function(x, y, plotit=TRUE, hollow=TRUE, xlab="X", ylab="Y", 
           trend="positive", clx=1, cly=1, 
           xthreshold=0.10, ythreshold=0.10,
           progress=TRUE, verbose=TRUE, listout=FALSE)
  {
  
  n = length(x)
  l = n-1
      dataset = data.frame(x=x, y=y,
                           xgroup=as.factor(c("a",rep("b", l))),
                           ygroup=as.factor(c("c",rep("d", l))),
                           Critical.x.value=as.numeric(rep(0.00, n)),
                           Sum.of.squares=as.numeric(rep(0.00, n)),
                           Critical.y.value=as.numeric(rep(0.00, n)),
                           q.i=as.integer(rep(0,n)),
                           q.ii=as.integer(rep(0,n)),
                           q.iii=as.integer(rep(0,n)),
                           q.iv=as.integer(rep(0,n)),
                           Q.i=as.integer(rep(0,n)),
                           Q.ii=as.integer(rep(0,n)),
                           Q.iii=as.integer(rep(0,n)),
                           Q.iv=as.integer(rep(0,n)),
                           Q.model=as.integer(rep(0,n)),
                           Q.err=as.integer(rep(0,n)),
                           Cramer.V=as.numeric(rep(0.00, n))
                           )
      dataset = dataset[with(dataset, order(x, y)),]
      
      for(i in (2:n)){
          dataset$Critical.x.value[i] = (dataset$x[i] + dataset$x[i-1])/2
          }
      dataset$Critical.x.value[1] = dataset$Critical.x.value[1]/2 

  m = n-2
  if(progress){tick=0}
  for(j in (3:m)){
    if(progress){tick=tick+1}
    if(progress){if(tick>(m-3)/100){tick=0;cat(".")}}
      for (i in (1:n)){
         dataset$xgroup[i] = "b"
         if(dataset$x[i] < dataset$Critical.x.value[j])
            dataset$xgroup[i] = as.factor("a")
      }
        fit = lm(y ~ xgroup, data=dataset)
        fit1 = anova(fit)
        dataset$Sum.of.squares[j] = (fit1[1,2])
      }
    if(progress){cat("\n")}
    
        dataset$Sum.of.squares[1] = min(dataset$Sum.of.squares[3:(n-2)])
        dataset$Sum.of.squares[2] = min(dataset$Sum.of.squares[3:(n-2)])
        dataset$Sum.of.squares[n-1] = min(dataset$Sum.of.squares[3:(n-2)])
        dataset$Sum.of.squares[n] = min(dataset$Sum.of.squares[3:(n-2)])
        max.ss = min(dataset$Sum.of.squares)+
                 (max(dataset$Sum.of.squares)-min(dataset$Sum.of.squares))*
                 (1-xthreshold)
        
              
      if(plotit){dataset4 = dataset[2:length(dataset[[1]]),]}

        dataset2 = dataset[with(dataset, order(-Sum.of.squares)),]
        rownames(dataset2) <- 1:nrow(dataset2)
        CLX = dataset2$Critical.x.value[clx]
            for (i in (1:n)){
               dataset$xgroup[i] = if(dataset$x[i] < CLX)  "a" else "b"
               }
  
  if(verbose){cat("\n") 
  cat("Critical x that maximize sum of squares:","\n","\n") 
  print(dataset2[dataset2$Sum.of.squares >= max.ss,][c("Critical.x.value", "Sum.of.squares")])}
  
  dataset = dataset[with(dataset, order(y, x)),]
  dataset$Critical.y.value[1] = 0
  for(k in c(2:n)){
     dataset$Critical.y.value[k] = (dataset$y[k]+dataset$y[k-1])/2
     }

  dataset$Critical.y.value[1] = min(dataset$Critical.y.value[2:n])

  if(progress){tick=0}
  for(j in c(1:n)){
    if(progress){tick=tick+1}
    if(progress){if(tick>(m-3)/100){tick=0;cat(".")}}
     for (i in c(1:n)){
       dataset$ygroup[i] = "d"
       if(dataset$y[i] < dataset$Critical.y.value[j])
          dataset$ygroup[i] = "c"
      }
      for (i in c(1:n)){
        dataset$q.i[i] = with(dataset, ifelse
        ((dataset$xgroup[i]=="a" & dataset$ygroup[i]=="d"), 1, 0))
        dataset$q.ii[i] = with(dataset, ifelse
        ((dataset$xgroup[i]=="b" & dataset$ygroup[i]=="d"), 1, 0))
        dataset$q.iii[i] = with(dataset, ifelse
        ((dataset$xgroup[i]=="b" & dataset$ygroup[i]=="c"), 1, 0))
        dataset$q.iv[i] = with(dataset, ifelse
        ((dataset$xgroup[i]=="a" & dataset$ygroup[i]=="c"), 1, 0)) 
       }
       dataset$Q.i[j]   = sum(dataset$q.i)
       dataset$Q.ii[j]  = sum(dataset$q.ii)
       dataset$Q.iii[j] = sum(dataset$q.iii)
       dataset$Q.iv[j]  = sum(dataset$q.iv)
       dataset$Cramer.V[j] = as.numeric(cramerV(matrix(c(dataset$Q.i[j],
                                                         dataset$Q.ii[j],
                                                         dataset$Q.iv[j],
                                                         dataset$Q.iii[j]),
                                                         nrow=2,
                                                         byrow=TRUE)))
       if (trend=="positive")
       dataset$Q.model[j] = dataset$Q.ii[j] + dataset$Q.iv[j]
       if (trend=="negative")
       dataset$Q.model[j] = dataset$Q.i[j] + dataset$Q.iii[j] 
       if (trend=="positive")
       dataset$Q.err[j] = dataset$Q.i[j] + dataset$Q.iii[j]
       if (trend=="negative")
       dataset$Q.err[j] = dataset$Q.ii[j] + dataset$Q.iv[j] 
     }
      if(progress){cat("\n")}
     
min.qerr = min(dataset$Q.err)+
          (max(dataset$Q.err)-min(dataset$Q.err))*(ythreshold)               

dataset3 = dataset[with(dataset, order(Q.err)),]

rownames(dataset3) <- 1:nrow(dataset3)
  if(verbose){cat("\n")
  cat("Critical y that minimize errors:","\n","\n")
  print(dataset3[dataset3$Q.err <= min.qerr,][c("Critical.y.value",
                   "Q.i","Q.ii","Q.iii","Q.iv",
                   "Q.model","Q.err","Cramer.V")])}

CLY = dataset3$Critical.y.value[cly]

 N = 1
 Z = data.frame(n  =as.numeric(rep(0.00, N)),
                CLx=as.numeric(rep(0.00, N)),
                SS =as.numeric(rep(0.00, N)),
                CLy=as.numeric(rep(0.00, N)),
                 Q.I=as.integer(rep(0,N)),
                 Q.II=as.integer(rep(0,N)),
                 Q.III=as.integer(rep(0,N)),
                 Q.IV=as.integer(rep(0,N)),
                 Q.Model=as.integer(rep(0.00,N)),
                 p.Model=as.numeric(rep(0.00,N)),
                 Q.Error=as.integer(rep(0.00,N)),
                 p.Error=as.numeric(rep(0.00,N)),
                 Fisher.p.value=as.numeric(rep(0.00,N)),
                 Cramer.V=as.numeric(rep(0.00,N))
                 )

  Z$n[1] = n
  Z$CLx[1] = CLX
  Z$SS[1] = dataset2$Sum.of.squares[1]
  Z$CLy[1] = CLY
  Z$Q.I[1] = dataset3$Q.i[cly]
  Z$Q.II[1] = dataset3$Q.ii[cly]
  Z$Q.III[1] = dataset3$Q.iii[cly]
  Z$Q.IV[1] = dataset3$Q.iv[cly]
  Z$Q.Model[1] = dataset3$Q.model[cly]
  Z$p.Model[1] = dataset3$Q.model[cly]/n
  Z$Q.Error[1] = dataset3$Q.err[cly]
  Z$p.Error[1] = dataset3$Q.err[cly]/n
  
  M=matrix(c(Z$Q.I[1],Z$Q.II[1],Z$Q.IV[1],Z$Q.III[1]),
           nrow=2,
           byrow=TRUE)
           
  Z$Fisher.p.value = fisher.test(M)$p.value
  Z$Cramer.V = dataset3$Cramer.V[cly]
    
    if (plotit)
     {
      if (trend=="positive")
       {
        pchi=as.integer(rep(16,n))
        pchi[(dataset$x<CLX)&(dataset$y<CLY)] = 16
        pchi[(dataset$x>CLX)&(dataset$y>CLY)] = 16
        pchi[(dataset$x>CLX)&(dataset$y<CLY)] = 1
        pchi[(dataset$x<CLX)&(dataset$y>CLY)] = 1
        }
     if (trend=="negative")
      {
       pchi=as.integer(rep(16,n))
       pchi[(dataset$x<CLX)&(dataset$y<CLY)] = 1
       pchi[(dataset$x>CLX)&(dataset$y>CLY)] = 1
       pchi[(dataset$x>CLX)&(dataset$y<CLY)] = 16
       pchi[(dataset$x<CLX)&(dataset$y>CLY)] = 16
       }
     if (hollow==FALSE)
      {
       pchi=as.integer(rep(16,n))
       }
    }         
    if (plotit) 
    plot(dataset$x,dataset$y,pch=16,xlab=xlab,ylab=ylab)           
    if (plotit){
    plot(Sum.of.squares~Critical.x.value,
         data=dataset4,
         xlab="Critical-x value",
         ylab="Sum of squares")}
    if (plotit) 
    plot(Q.err~Critical.y.value, data=dataset, xlab="Critical-y value",
         ylab="Number of points in error quadrants")
    if (plotit) 
    plot(Cramer.V~Critical.y.value, data=dataset, xlab="Critical-y value",
         ylab="Cramer's V")
    if (plotit)
    {
     plot(dataset$x,dataset$y,pch=pchi,
          xlab=xlab,
          ylab=ylab)
     abline(v=CLX, col="blue")
     abline(h=CLY, col="blue")
     max.x = max(dataset$x)
     max.y = max(dataset$y)
     min.x = min(dataset$x)
     min.y = min(dataset$y)
     text (min.x+(max.x-min.x)*0.01, min.y+(max.y-min.y)*0.99,
           labels="I")
     text (min.x+(max.x-min.x)*0.99, min.y+(max.y-min.y)*0.99,
           labels="II")
     text (min.x+(max.x-min.x)*0.99, min.y+(max.y-min.y)*0.01,
           labels="III")
     text (min.x+(max.x-min.x)*0.01, min.y+(max.y-min.y)*0.01,
           labels="IV")
     }
 
if(verbose){cat("\n")    
cat("n         = Number of observations","\n")
cat("CLx       = Critical value of x","\n")
cat("SS        = Sum of squares for that critical value of x","\n")
cat("CLy       = Critical value of y","\n")
cat("Q         = Number of observations which fall into quadrants I, II, III, IV",
          "\n")
cat("Q.Model   = Total observations which fall into the quadrants predicted by the model",
          "\n")
cat("p.Model   = Percent observations which fall into the quadrants predicted by the model",
          "\n")
cat("Q.Error   = Observations which do not fall into the quadrants predicted by the model",
          "\n")
cat("p.Error   = Percent observations which do not fall into the quadrants predicted by the model",
          "\n")
cat("Fisher.p  = p-value from Fisher exact test dividing data into these quadrants",
          "\n")
cat("Cramer.V  = Cramer's V statistic from dividing data into these quadrants",
          "\n")
cat("\n") 
cat("Final model:","\n","\n")}
if(!listout){
  return(Z)}

  if(listout){
  print(Z)
  cat("\n") 
  Zed = list(Critical.x = dataset2[c("Critical.x.value", "Sum.of.squares")],
             Critical.y = dataset3[c("Critical.y.value",
                                     "Q.i","Q.ii","Q.iii","Q.iv",
                                     "Q.model","Q.err","Cramer.V")],
             Final.model = Z)
  return(Zed)}
}