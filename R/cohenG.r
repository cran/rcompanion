#' @title Cohen's g and odds ratio for paired contingency tables
#'
#' @description Calculates Cohen's g and odds ratio
#'              for paired contingency tables, such as those that 
#'              might be analyzed with   
#'              McNemar or McNemar-Bowker tests.
#' 
#' @param x A two-way contingency table. It must be square. 
#'          It can have two or
#'          more levels for each dimension.
#' @param ci If \code{TRUE}, returns confidence intervals by bootstrap.
#'           May be slow.
#' @param conf The level for the confidence interval.
#' @param type The type of confidence interval to use.
#'             Can be any of "\code{norm}", "\code{basic}", 
#'                           "\code{perc}", or "\code{bca}".
#'             Passed to \code{boot.ci}.
#' @param R The number of replications to use for bootstrap.
#' @param histogram If \code{TRUE}, produces a histogram of bootstrapped values.
#' @param digits The number of significant digits in the output.
#' @param reportIncomplete If \code{FALSE} (the default),
#'                         \code{NA} will be reported in cases where there
#'                         are instances of the calculation of the statistic
#'                         failing during the bootstrap procedure. 
#' @param ... Additional arguments (ignored).
#' 
#' @details For a 2 x 2 table, where a and d are the concordant cells
#'          and b and c are discordant cells:
#'          Odds ratio is b/c;
#'          P is b/(b+c);
#'          and Cohen's g is P - 0.5.
#'          
#'          In the 2 x 2 case, the statistics are directional.
#'          That is, when cell [1, 2] in the table is greater than
#'          cell [2, 1], OR is greater than 1, P is greater than 0.5,
#'          and g is positive.  
#'          
#'          In the opposite case, OR is less than 1,
#'          P is less than 0.5, and g is negative.
#'          
#'          In the 2 x 2 case, when the effect is small, the
#'          confidence interval for OR can pass through 1,
#'          for g can pass through 0, and for P can pass through 0.5.
#'          
#'          For tables larger than 2 x 2, the statistics are not directional.
#'          That is, OR is always >= 1, P is always >= 0.5, and
#'          g is always positive.  
#'          Because of this, if \code{type="perc"}, the confidence interval will
#'          never cross the values for no effect 
#'          (OR = 1, P = 0.5, or g = 0).
#'          Because of this, the confidence interval range
#'          in this case should not
#'          be used for statistical inference.
#'          However, if \code{type="norm"}, the confidence interval
#'          may cross the values for no effect.
#'          
#'           When the reported statistics are close to their extremes,
#'           or with small counts, 
#'           the confidence intervals 
#'           determined by this
#'           method may not be reliable, or the procedure may fail.
#'        
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/H_05.html}
#' @concept Effect size Cohen g McNemar Bowker contingency symmetry nominal
#' @return A list containing: a data frame of results of the global statistics;
#'         and a data frame of results of the pairwise statistics.
#'         
#' @seealso \code{\link{nominalSymmetryTest}}
#'         
#' @examples
#' ### 2 x 2 repeated matrix example
#' data(AndersonRainBarrel)
#' cohenG(AndersonRainBarrel)
#'                     
#' ### 3 x 3 repeated matrix
#' data(AndersonRainGarden)
#' cohenG(AndersonRainGarden)
#' 
#' @export

cohenG = 
  function(x, ci=FALSE, conf=0.95, type="perc",
           R=1000, histogram=FALSE, 
           digits=3, reportIncomplete=FALSE, ...) {
    
     n = nrow(x)
     m = ncol(x)
     if((n < 2) | (m != n)){
        stop("Matrix must be square with at least two rows and columns")}
     N = n*n-n
     Y = rep(1/N, N)
     X = rep(NA, N)
     k=0
     for(i in 1:(n)){
        for(j in 1:n){
        if(i != j){
           k=k+1
           X[k]=x[i,j]
          }
        }
     }
    if(n==2){
      b  = x[1,2]
      c  = x[2,1]
      OR = b/c
      OR = signif(OR, digits=digits)
      P  = b/(b+c)
      P  = signif(P, digits=digits)
      g  = P - 0.5
      g  = signif(g, digits=digits)
    }

  Dimensions = paste(n, "x", n)
  
  if(n==2 & ci==FALSE){
     V = data.frame(Dimensions, OR, P, g)
     W = list(V)
     names(W) = c("Global.statistics")
     return(W)
  }
  
    if(n==2 & ci==TRUE){
     Statistic   = c("OR", "P", "g")
     Value       = c(OR, P, g)
     Dimensions  = rep(Dimensions, 3)
     
          if(is.nan(OR) | is.nan(P) | is.nan(g)){
              Value    = c(NaN, NaN, NaN)
              lower.ci = c(NA, NA, NA)
              upper.ci = c(NA, NA, NA)
              V        = data.frame(Dimensions, Statistic, Value, 
                                   lower.ci, upper.ci)
              W        = list(V)
              names(W) = c("Global.statistics")
              return(W)
              }
     
    X = as.table(x)
    Counts = as.data.frame(X)
    Long = Counts[rep(row.names(Counts), Counts$Freq), c(1, 2)]
    rownames(Long) = seq(1:nrow(Long))
    
    L1     = length(unique(droplevels(Long[,1])))
    L2     = length(unique(droplevels(Long[,2])))
    
    Function = function(input, index){
        Input = input[index,]
                
        NOTEQUAL=0
        if(length(unique(droplevels(Input[,1]))) != L1 |
           length(unique(droplevels(Input[,2]))) != L2){NOTEQUAL=1}
             
        if(NOTEQUAL==1){FLAG=1; return(c(NA,FLAG))}
             
        if(NOTEQUAL==0){
                
                Table = table(Input)
                b  = Table[1,2]
                c  = Table[2,1]
                Stat = b/c
                FLAG = 0
                return(c(Stat,FLAG))}
    }
    
    Boot = boot(Long, Function, R=R)
    BCI  = boot.ci(Boot, conf=conf, type=type)
    if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3]}
    if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5]}
    if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5]}
    if(type=="bca")  {CI1=BCI$bca[4];     CI2=BCI$bca[5]}
    
    if(sum(Boot$t[,2])>0 & reportIncomplete==FALSE) {CI1=NA; CI2=NA}
  
    OR1=signif(CI1, digits=digits)
    OR2=signif(CI2, digits=digits)
    
    if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray", xlab="OR", main="")}
    
    Function = function(input, index){
        Input = input[index,]
                
        NOTEQUAL=0
        if(length(unique(droplevels(Input[,1]))) != L1 |
           length(unique(droplevels(Input[,2]))) != L2){NOTEQUAL=1}
             
        if(NOTEQUAL==1){FLAG=1; return(c(NA,FLAG))}
             
        if(NOTEQUAL==0){
                Table = table(Input)
                b  = Table[1,2]
                c  = Table[2,1]
                Stat = b/(b+c)
                FLAG = 0
                return(c(Stat,FLAG))}
    }
    
    Boot = boot(Long, Function, R=R)
    BCI  = boot.ci(Boot, conf=conf, type=type)
    if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3]}
    if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5]}
    if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5]}
    if(type=="bca")  {CI1=BCI$bca[4];     CI2=BCI$bca[5]}
    
    if(sum(Boot$t[,2])>0 & reportIncomplete==FALSE) {CI1=NA; CI2=NA}
    
    P1=signif(CI1, digits=digits)
    P2=signif(CI2, digits=digits)
    
    g1 = P1-0.5
    g2 = P2-0.5
  
    if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray", xlab="P", main="")}
    if(histogram==TRUE){hist(Boot$t[,1]-0.5, col = "darkgray", 
                                             xlab="g", main="")}
    
     lower.ci    = c(OR1, P1, g1)
     upper.ci    = c(OR2, P2, g2)
     
     V = data.frame(Dimensions, Statistic, Value, lower.ci, upper.ci)
     W = list(V)
     names(W) = c("Global.statistics")
     return(W)
  }
  
  if(n>2){
     N = n*(n-1)/2
     Z = data.frame(Comparison = rep("A", N),
                    OR         = rep(NA, N),
                    P          = rep(NA, N),
                    g          = rep(NA, N),
                    stringsAsFactors=FALSE)
     k=0
     NUMERATOR = 0
     DENOMINATOR = 0
     DENOMINATOR1 = 0
     for(i in 1:(n-1)){
        for(j in (i+1):n){
           k=k+1
           Namea = as.character(rownames(x)[i])
           Nameb = as.character(colnames(x)[i])
           Namec = as.character(rownames(x)[j])
           Named = as.character(colnames(x)[j])
           b = x[i,j]
           c = x[j,i]
           OR = max(b/c, c/b)
           P  = max(b/(b+c), c/(c+b))
           g  = P - 0.5
           OR = signif(OR, digits=digits)
           P  = signif(P, digits=digits)
           g  = signif(g, digits=digits)
           Z[k,] =c( paste0(Namea, "/", Nameb, " : ", Namec, "/", Named),
                    OR, P, g)
           NUMERATOR   = NUMERATOR + max(b, c)
           DENOMINATOR = DENOMINATOR + min(b, c)
           DENOMINATOR1 = DENOMINATOR1 + b + c
        }
     }
     OR = NUMERATOR / DENOMINATOR
     P  = NUMERATOR / DENOMINATOR1
     g  = P - 0.5
     OR = signif(OR, digits=digits)
     P  = signif(P, digits=digits)
     g  = signif(g, digits=digits)
     
     V  = data.frame(Dimensions, OR, P, g)
  }
     
  if(n>2 & ci==TRUE){
    
           if(is.nan(OR) | is.nan(P) | is.nan(g)){
               Statistic   = c("OR", "P", "g")
               Value    = c(NaN, NaN, NaN)
               lower.ci = c(NA, NA, NA)
               upper.ci = c(NA, NA, NA)
               V        = data.frame(Dimensions, Statistic, Value, 
                                     lower.ci, upper.ci)
               W        = list(V, Z)
               names(W) = c("Global.statistics", "Pairwise.statistics")
               return(W)
              }

         Funny = function(x){
           k=0
           NUMERATOR = 0
           DENOMINATOR = 0
           DENOMINATOR1 = 0
           for(i in 1:(n-1)){
            for(j in (i+1):n){
             k=k+1
             b = x[i,j]
             c = x[j,i]
             NUMERATOR   = NUMERATOR + max(b, c)
             DENOMINATOR = DENOMINATOR + min(b, c)
             DENOMINATOR1 = DENOMINATOR1 + b + c
             }
           }
          OR = NUMERATOR / DENOMINATOR
          P  = NUMERATOR / DENOMINATOR1
          g  = P - 0.5
          Output = c(OR, P, g)
          names(Output) = c("OR", "P", "g")
          return(Output)
         }
  
  X = as.table(x)
    Counts = as.data.frame(X)
    Long = Counts[rep(row.names(Counts), Counts$Freq), c(1, 2)]
    rownames(Long) = seq(1:nrow(Long))
    
        L1     = length(unique(droplevels(Long[,1])))
        L2     = length(unique(droplevels(Long[,2])))
    
    Function = function(input, index){
                Input = input[index,]
                
          NOTEQUAL=0
          if(length(unique(droplevels(Input[,1]))) != L1 |
             length(unique(droplevels(Input[,2]))) != L2){NOTEQUAL=1}
             
          if(NOTEQUAL==1){FLAG=1; return(c(NA,FLAG))}
             
          if(NOTEQUAL==0){
                Table = table(Input)
                Stat = Funny(Table)[1]
                FLAG=0
                return(c(Stat,FLAG))}
    }
    
    Boot = boot(Long, Function, R=R)
    BCI  = boot.ci(Boot, conf=conf, type=type)
    if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3]}
    if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5]}
    if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5]}
    if(type=="bca")  {CI1=BCI$bca[4];     CI2=BCI$bca[5]}  
  
    if(sum(Boot$t[,2])>0 & reportIncomplete==FALSE) {CI1=NA; CI2=NA}
    
    OR1=signif(CI1, digits=digits)
    OR2=signif(CI2, digits=digits)
    
    if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray", xlab="OR", main="")}
    
  Function = function(input, index){
        Input = input[index,]
                
        NOTEQUAL=0
        if(length(unique(droplevels(Input[,1]))) != L1 |
           length(unique(droplevels(Input[,2]))) != L2){NOTEQUAL=1}
        
        if(NOTEQUAL==1){FLAG=1; return(c(NA,FLAG))}
                
        if(NOTEQUAL==0){
                Table = table(Input)
                Stat = Funny(Table)[2]
                FLAG = 0
                return(c(Stat,FLAG))}
  }
    
    Boot = boot(Long, Function, R=R)
    BCI  = boot.ci(Boot, conf=conf, type=type)
    if(type=="norm") {CI1=BCI$normal[2];  CI2=BCI$normal[3]}
    if(type=="basic"){CI1=BCI$basic[4];   CI2=BCI$basic[5]}
    if(type=="perc") {CI1=BCI$percent[4]; CI2=BCI$percent[5]}
    if(type=="bca")  {CI1=BCI$bca[4];     CI2=BCI$bca[5]}
    
    if(sum(Boot$t[,2])>0 & reportIncomplete==FALSE) {CI1=NA; CI2=NA}
  
    P1=signif(CI1, digits=digits)
    P2=signif(CI2, digits=digits)
    
    g1 = P1-0.5
    g2 = P2-0.5
  
    if(histogram==TRUE){hist(Boot$t[,1], col = "darkgray", xlab="P", main="")}
    if(histogram==TRUE){hist(Boot$t[,1]-0.5, col = "darkgray",
                                             xlab="g", main="")}
    
    Statistic   = c("OR", "P", "g")
    Value       = c(OR, P, g)
    Dimensions  = rep(Dimensions, 3)
    lower.ci    = c(OR1, P1, g1)
    upper.ci    = c(OR2, P2, g2)
    V = data.frame(Dimensions, Statistic, Value, lower.ci, upper.ci)
  }

     W = list(V, Z)
     names(W) = c("Global.statistics",
                  "Pairwise.statistics")
     return(W)
  }
