#' @title Correlation and measures of association
#'
#' @description Produces measures of association for all variables
#'               in a data frame with confidence intervals when available.
#' 
#' @param data         A data frame.
#' 
#' @param ci           If \code{TRUE}, calculates confidence intervals 
#'                     for methods requiring bootstrap.  
#'                     If \code{FALSE}, will return only those
#'                     confidence intervals from methods not 
#'                     requiring bootstrap.
#'                     
#' @param conf         The confidence level for confidence intervals.
#' 
#' @param printClasses If \code{TRUE}, prints a table of classes for
#'                     all variables.
#'                     
#' @param progress     If \code{TRUE}, prints progress bar when bootstrap
#'                     methods are called.
#'                     
#' @param methodNum    The method for the correlation for two numeric variables.
#'                     The default is \code{"pearson"}. Other options are
#'                     \code{"spearman"} and \code{"kendall"}.
#'                     
#' @param methodOrd    The method for the correlation for two ordinal variables.
#'                     The default is \code{"kendall"}, with Kendall's tau-c
#'                     used. Other option is \code{"spearman"}.
#'                     
#' @param methodNumOrd The method for the correlation of a numeric and
#'                     an ordinal variable.
#'                     The default is \code{"pearson"}. Other options are
#'                     \code{"spearman"} and \code{"kendall"}.
#'                     
#' @param methodNumNom The method for the correlation of a numeric and
#'                     a nominal variable.
#'                     
#'                     The default is \code{"eta"}, which is the square root
#'                     of the r-squared value from anova. 
#'                     The other option is \code{"epsilon"}, which is the same,
#'                     except with the numeric value rank-transformed.
#'                     
#' @param methodNumBin The method for the correlation of a numeric and
#'                     a binary variable.
#'                     The default is \code{"pearson"}.
#'                     The other option is \code{"glass"}, which uses the 
#'                     Glass rank biserial correlation.
#'                     
#' @param testChisq    The method for the test of two nominal variables.
#'                     The default is \code{"chisq"}. The other option is
#'                     \code{"fisher"}.
#'                     
#' @param R            The number of replications to use for bootstrap 
#'                     confidence intervals for applicable methods.
#'                     
#' @param reportIncomplete If \code{FALSE}, \code{NA} will be reported in cases 
#'                         where there are instances of the calculation of the
#'                         statistic failing during the bootstrap procedure.
#'                         
#' @param  na.action   If \code{"na.omit"}, the function will use only 
#'                     complete cases, assessed on a bivariate basis.  
#'                     The other option is \code{"na.pass"}.
#'             
#' @param correct      Passed to \code{chisq.test}.
#' 
#' @param  digits     The number of decimal places in the output of most 
#'                    statistics.
#'                     
#' @param  pDigits   The number of decimal places in the output for 
#'                    p-values.          
#'              
#' @param ... Other arguments.
#' 
#' @details  It’s important that variables are assigned the correct class
#'           to get an appropriate measure of association.  
#'           That is, factor variables should be of class "factor", 
#'           not "character". Ordered factors should be ordered factors
#'           (and have their levels in the correct order!).
#'             
#'           Date variables are treated as numeric.
#'           
#'           The default for measures of association tend to be
#'           "parametric" type. That is, e.g. Pearson correlation
#'           where appropriate.
#'           
#'           Nonparametric measures of association will be reported
#'           with the options
#'           \code{methodNum = "spearman", methodNumNom = "epsilon",
#'           methodNumBin = "glass", methodNumOrd="spearman"}.
#'
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{https://rcompanion.org/handbook/I_14.html}
#' 
#' @seealso \code{\link{phi}}, 
#'          \code{\link{spearmanRho}},
#'          \code{\link{cramerV}},
#'          \code{\link{freemanTheta}},
#'          \code{\link{wilcoxonRG}}
#'
#' @concept correlation
#' @concept association
#' @concept Pearson
#' @concept Kendall
#' @concept Spearman
#' @concept phi
#' @concept Cramer's V
#' @concept Freemans's theta
#' @concept Glass rank biserial coefficient
#' @concept eta
#' @concept epsilon
#' 
#' @return A data frame of variables, association statistics, p-values, and
#'          confidence intervals.
#'         
#' @examples
#' 
#' Length   = c(0.29, 0.25, NA, 0.40, 0.50, 0.57, 0.62, 0.88, 0.99, 0.90)
#' Rating   = factor(ordered=TRUE, levels=c("Low", "Medium", "High"),
#'                   x = rep(c("Low", "Medium", "High"), c(3,3,4)))
#' Color    = factor(rep(c("Red", "Green", "Blue"), c(4,4,2)))
#' Flag     = factor(rep(c(TRUE, FALSE, TRUE), c(5,4,1)))
#' Answer   = factor(rep(c("Yes", "No", "Yes"), c(4,3,3)), levels=c("Yes", "No"))
#' Location = factor(rep(c("Home", "Away", "Other"), c(2,4,4)))
#' Distance = factor(ordered=TRUE, levels=c("Low", "Medium", "High"),
#'                   x = rep(c("Low", "Medium", "High"), c(5,2,3))) 
#' Start    = seq(as.Date("2024-01-01"), by = "month", length.out = 10)
#' Data = data.frame(Length, Rating, Color, Flag, Answer, Location, Distance, Start)  
#' correlation(Data)
#' 
#' 
#' 
#' @importFrom stats      na.omit cor.test chisq.test fisher.test lm
#' @importFrom rcompanion spearmanRho phi cramerV freemanTheta wilcoxonRG
#' @importFrom coin       lbl_test chisq_test pvalue
#' @importFrom DescTools  StuartTauC KendallTauB
#' 
#' @export

correlation = 
function (data=NULL, 
          printClasses = FALSE,
          progress     = TRUE,
          methodNum    = "pearson",
          methodOrd    = "kendall",
          methodNumOrd = "spearman",
          methodNumNom = "eta",
          methodNumBin = "pearson",
          testChisq    = "chisq",
          ci           = FALSE,
          conf         = 0.95,
          R            = 1000,
          correct      = FALSE,
          reportIncomplete = TRUE,
          na.action    = "na.omit",
          digits       = 3,
          pDigits      = 4,
          ...) {
 
 m = ncol(data)
 
 Y = data.frame(Variable  = as.character(rep(NA, m)),
                Class     = as.character(rep(NA, m)),
                Treatment = as.character(rep(NA, m)))
 
  for(l in 1:m){
   
   Y[l,1] = colnames(data)[l]
   Y[l,2] = class(data[[l]])[1]
   
   Class = Y[l,2]
   
   Treatment = "Character"
  
   if(Class == "integer"){Treatment="Numeric"}
   if(Class == "numeric"){Treatment="Numeric"}
   if(Class == "ordered"){Treatment="Ordinal"}
   if(Class == "factor") {Treatment="Nominal"}
   if(Class == "logical"){Treatment="Binary"}
   
   if(Class == "factor"){data[,l] = droplevels(data[,l])}
   if(Class == "factor" & length(levels(data[,l]))==2){Treatment="Binary"}

   if(Class == "Date"){
    Treatment="Numeric"
    data[,l] = as.numeric(data[,l])
    }
      
   Y[l,3] = Treatment
  } # End m loop
 
 if(printClasses){
  cat("\n")
  print(Y)
  cat("\n")
  }

 n = m*(m-1)/2
 
 Z = data.frame(Var1      = as.character(rep(NA, n)),
                Var2      = as.character(rep(NA, n)),
                Type      = as.character(rep(NA, n)),
                N         = as.numeric  (rep(NA, n)),
                Measure   = as.character(rep(NA, n)),
                Statistic = as.numeric  (rep(NA, n)),
                Lower.CL  = as.numeric  (rep(NA, n)),
                Upper.CL  = as.numeric  (rep(NA, n)),
                Test      = as.character(rep(NA, n)),
                p.value   = as.numeric  (rep(NA, n)),
                Signif    = as.character(rep(NA, n))
                )
k=0               
 for(i in 1:(m-1)){
  for(j in (i+1):m){
   k=k+1
   
   if(ci){if(progress){cat(".")}}
   
   A = data[,i]
   B = data[,j]
            
   if(na.action=="na.omit"){
    NAH = na.omit(data.frame(A,B))
    A = NAH$A
    B = NAH$B
   }
   
   Z[k,1] = colnames(data)[i]
   Z[k,2] = colnames(data)[j]
   Z[k,3] = paste0(Y[i,3]," x ",Y[j,3])
   Z[k,4] = length(A)
   
   ### START Numeric and Numeric
   
   if( Y[i,3] == "Numeric" & Y[j,3] == "Numeric" ){
      if(methodNum == "pearson"){
        Z[k,5] = "Pearson"
        Test   = suppressWarnings(cor.test(~ A + B, method="pearson", conf.level=conf))
        Z[k,6] = round(Test$estimate, digits)
        Z[k,7] = round(Test$conf.int[1], digits)
        Z[k,8] = round(Test$conf.int[2], digits)
        Z[k,9] = "cor.test"
        Z[k,10] = round(Test$p.value, pDigits)
      }
      if(methodNum == "spearman"){
       Z[k,5] = "Spearman"
       Test   = suppressWarnings(cor.test(~ rank(A) + rank(B), method="pearson", conf.level=conf))
       Z[k,6] = round(Test$estimate, digits)
       Z[k,7] = round(Test$conf.int[1], digits)
       Z[k,8] = round(Test$conf.int[2], digits)
       Z[k,9] = "cor.test"
       Z[k,10] = round(Test$p.value, pDigits)
      }
     if(methodNum == "kendall"){
       Z[k,5] = "Kendall"
       Test   = suppressWarnings(spearmanRho(x=A, y=B, method="kendall", ci=ci, R=R, reportIncomplete=reportIncomplete, conf=conf))
       if(!ci){Z[k,6] = round(Test, digits)}
       if(ci) {Z[k,6] = round(Test$tau, digits)}
       if(ci) {Z[k,7] = round(Test$lower.ci, digits)}
       if(ci) {Z[k,8] = round(Test$upper.ci, digits)}
       Test2  = suppressWarnings(cor.test(~ A + B, method="kendall"))
       Z[k,9] = "cor.test"
       Z[k,10] = round(Test2$p.value, pDigits)
      }
   } # End numeric x numeric
   
  ### START Ordinal and Ordinal
   
  if( Y[i,3] == "Ordinal" & Y[j,3] == "Ordinal" ){
    if(methodOrd == "spearman"){
     Z[k,5] = "Spearman"
     Test   = suppressWarnings(cor.test(~ rank(as.numeric(A)) + rank(as.numeric(B)), method="pearson", conf.level=conf))
     Z[k,6] = round(Test$estimate, digits)
     Z[k,7] = round(Test$conf.int[1], digits)
     Z[k,8] = round(Test$conf.int[2], digits)
     Z[k,9] = "cor.test"
     Z[k,10] = round(Test$p.value, pDigits)
    }
    if(methodOrd == "kendall"){
       Z[k,5] = "Kendall"
       Tab = table(A,B)
       Test   = suppressWarnings(StuartTauC(Tab, conf.level=conf))
       Z[k,6] = round(Test[1], digits)
       Z[k,7] = round(Test[2], digits)
       Z[k,8] = round(Test[3], digits)
       Test2 = suppressWarnings(lbl_test(Tab))
       Z[k,9] = "Linear by linear"
       Z[k,10] = round(pvalue(Test2), pDigits)
    }
   } # End ordinal x ordinal
   
   ### START Binary and Binary
   
   if( Y[i,3] == "Binary" & Y[j,3] == "Binary" ){
    Z[k,5] = "Phi"
    Tab = table(A,B)
    Test   = suppressWarnings(phi(Tab, ci=ci, R=R, reportIncomplete=reportIncomplete, conf=conf))
    if(!ci){Z[k,6] = round(Test, digits)}
    if(ci) {Z[k,6] = round(Test$phi, digits)}
    if(ci) {Z[k,7] = round(Test$lower.ci, digits)}
    if(ci) {Z[k,8] = round(Test$upper.ci, digits)}
    if(testChisq=="chisq"){
         Test2 = suppressWarnings(chisq.test(Tab, correct=correct))
         Z[k,9] = "chisq.test"
         Z[k,10] = round(Test2$p.value, pDigits)
         }
    if(testChisq=="fisher"){
         Test2 = suppressWarnings(fisher.test(Tab))
         Z[k,9] = "fisher.test"
         Z[k,10] = round(Test2$p.value, pDigits)
        }
       } # End binary and binary
   
   ### START Nominal and Nominal
   
   if( Y[i,3] == "Nominal" & Y[j,3] == "Nominal" |
       Y[i,3] == "Nominal" & Y[j,3] == "Binary"  |
       Y[i,3] == "Binary"  & Y[j,3] == "Nominal"){
      Z[k,5] = "Cramer"
      Tab = table(A,B)
      Test   = suppressWarnings(cramerV(Tab, ci=ci, R=R, reportIncomplete=reportIncomplete, conf=conf))
      if(!ci){Z[k,6] = round(Test, digits)}
      if(ci) {Z[k,6] = round(Test$Cramer.V, digits)}
      if(ci) {Z[k,7] = round(Test$lower.ci, digits)}
      if(ci) {Z[k,8] = round(Test$upper.ci, digits)}
     if(testChisq=="chisq"){
      Test2 = suppressWarnings(chisq.test(Tab, correct=correct))
       Z[k,9] = "chisq.test"
       Z[k,10] = round(Test2$p.value, pDigits)
      }
    if(testChisq=="fisher"){
     Test2 = suppressWarnings(fisher.test(Tab))
     Z[k,9] = "fisher.test"
     Z[k,10] = round(Test2$p.value, pDigits)
    }
   } # End nominal x nominal
   
   ### START Numeric and Ordinal
   
  if( Y[i,3] == "Ordinal" & Y[j,3] == "Numeric" ){
    if(methodNumOrd == "pearson"){
     Z[k,5] = "Pearson"
     Test   = suppressWarnings(cor.test(~ as.numeric(A) + B, method="pearson", conf.level=conf))
     Z[k,6] = round(Test$estimate, digits)
     Z[k,7] = round(Test$conf.int[1], digits)
     Z[k,8] = round(Test$conf.int[2], digits)
     Z[k,9] = "cor.test"
     Z[k,10] = round(Test$p.value, pDigits)
    }
    if(methodNumOrd == "spearman"){
     Z[k,5] = "Spearman"
     Test   = suppressWarnings(cor.test(~ rank(as.numeric(A)) + rank(B), method="pearson", conf.level=conf))
     Z[k,6] = round(Test$estimate, digits)
     Z[k,7] = round(Test$conf.int[1], digits)
     Z[k,8] = round(Test$conf.int[2], digits)
     Z[k,9] = "cor.test"
     Z[k,10] = round(Test$p.value, pDigits)
    }
    if(methodNumOrd == "kendall"){
     Z[k,5] = "Kendall"
     Test   = suppressWarnings(KendallTauB(x=A, y=B, conf.level=conf))
     Z[k,6] = round(Test[1], digits)
     Z[k,7] = round(Test[2], digits)
     Z[k,8] = round(Test[3], digits)
     Test2  = suppressWarnings(cor.test(~ as.numeric(A) + B, method="pearson"))
     Z[k,9] = "cor.test"
     Z[k,10] = round(Test2$p.value, pDigits)
     }
   }
   
   if( Y[i,3] == "Numeric" & Y[j,3] == "Ordinal" ){
    if(methodNumOrd == "pearson"){
     Z[k,5] = "Pearson"
     Test   = suppressWarnings(cor.test(~ A + as.numeric(B), method="pearson", conf.level=conf))
     Z[k,6] = round(Test$estimate, digits)
     Z[k,7] = round(Test$conf.int[1], digits)
     Z[k,8] = round(Test$conf.int[2], digits)
     Z[k,9] = "cor.test"
     Z[k,10] = round(Test$p.value, pDigits)
    }
    if(methodNumOrd == "spearman"){
     Z[k,5] = "Spearman"
     Test   = suppressWarnings(cor.test(~ rank(A) + rank(as.numeric(B)), method="pearson", conf.level=conf))
     Z[k,6] = round(Test$estimate, digits)
     Z[k,7] = round(Test$conf.int[1], digits)
     Z[k,8] = round(Test$conf.int[2], digits)
     Z[k,9] = "cor.test"
     Z[k,10] = round(Test$p.value, pDigits)
    }
    if(methodNumOrd == "kendall"){
     Z[k,5] = "Kendall"
     Test   = suppressWarnings(KendallTauB(x=A,y=B, conf.level=conf))
     Z[k,6] = round(Test[1], digits)
     Z[k,7] = round(Test[2], digits)
     Z[k,8] = round(Test[3], digits)
     Z[k,9] = "cor.test"
     Test2  = suppressWarnings(cor.test(~ A + as.numeric(B), method="pearson"))
     Z[k,10] = round(Test2$p.value, pDigits)
    }
   } # End ordinal x numeric
   
   ### START Nominal and Ordinal
   
   if( Y[i,3] == "Nominal" & Y[j,3] == "Ordinal" ){
    Z[k,5] = "Freeman"
    Tab = table(A, B)
    names(dimnames(Tab)) = c("Nominal", "Ordinal")
    Test   = suppressWarnings(freemanTheta(Tab, ci=ci, R=R, reportIncomplete=reportIncomplete, conf=conf))
    if(!ci){Z[k,6] = round(Test, digits)}
    if(ci) {Z[k,6] = round(Test$Freeman.theta, digits)}
    if(ci) {Z[k,7] = round(Test$lower.ci, digits)}
    if(ci) {Z[k,8] = round(Test$upper.ci, digits)}
    Test2 = suppressWarnings(chisq_test(Tab, scores = list("Ordinal" = 1:ncol(Tab))))
    Z[k,9] = "Cochran-Armitage"
    Z[k,10] = round(pvalue(Test2), pDigits)
   }
   
   if( Y[i,3] == "Ordinal" & Y[j,3] == "Nominal" ){
    Z[k,5] = "Freeman"
    Tab = table(B, A)
    names(dimnames(Tab)) = c("Nominal", "Ordinal")
    Test   = suppressWarnings(freemanTheta(Tab, ci=ci, R=R, reportIncomplete=reportIncomplete, conf=conf))
    if(!ci){Z[k,6] = round(Test, digits)}
    if(ci) {Z[k,6] = round(Test$Freeman.theta, digits)}
    if(ci) {Z[k,7] = round(Test$lower.ci, digits)}
    if(ci) {Z[k,8] = round(Test$upper.ci, digits)}
    Test2 = suppressWarnings(coin::chisq_test(Tab, scores = list("Ordinal" = 1:ncol(Tab))))
    Z[k,9] = "Cochran-Armitage"
    Z[k,10] = round(pvalue(Test2), pDigits)
   }
   # End nominal and ordinal
   
   ### START Numeric and Binary
   
   if( Y[i,3] == "Binary" & Y[j,3] == "Numeric" ){
    if(methodNumBin == "pearson"){
     Z[k,5] = "Pearson"
     Test = suppressWarnings(cor.test(~ as.numeric(A) + B, conf.level=conf))
     Z[k,6] = round(Test$estimate, digits)
     Z[k,7] = round(Test$conf.int[1], digits)
     Z[k,8] = round(Test$conf.int[2], digits)
     Z[k,9] = "cor.test"
     Z[k,10] = round(Test$p.value, pDigits)
    }
    if(methodNumBin == "glass"){
     Z[k,5] = "Glass rank biserial"
     Test = suppressWarnings(wilcoxonRG(x = B, g = A, ci=ci, R=R, reportIncomplete=reportIncomplete, conf=conf))
     if(!ci){
      Z[k,6] = round(-1 * Test, digits)
      }
     if(ci){
       Z[k,6] = round(-1 * Test$rg, digits)
       Z[k,7] = round(-1 * Test$lower.ci, digits)
       Z[k,8] = round(-1 * Test$upper.ci, digits)
     }
     Z[k,9] = "wilcox.test"
     Test2  = suppressWarnings(wilcox.test(B ~ A))
     Z[k,10] = round(Test2$p.value, pDigits)
    }
   }
   
   if( Y[i,3] == "Numeric" & Y[j,3] == "Binary" ){
    if(methodNumBin == "pearson"){
     Test = suppressWarnings(cor.test(~ A + as.numeric(B), conf.level=conf))
     Z[k,5] = "Pearson"
     Z[k,6] = round(Test$estimate, digits)
     Z[k,7] = round(Test$conf.int[1], digits)
     Z[k,8] = round(Test$conf.int[2], digits)
     Z[k,9] = "cor.test"
     Z[k,10] = round(Test$p.value, pDigits)
    }
    if(methodNumBin == "glass"){
     Z[k,5] = "Glass rank biserial"
     Test = suppressWarnings(wilcoxonRG(x = A, g = B, ci=ci, R=R, reportIncomplete=reportIncomplete, conf=conf))
     if(!ci){
      Z[k,6] = round(-1 * Test, digits)
      }
     if(ci){
      Z[k,6] = round(-1 * Test$rg, digits)
      Z[k,7] = round(-1 * Test$lower.ci, digits)
      Z[k,8] = round(-1 * Test$upper.ci, digits)
     }
     Z[k,9] = "wilcox.test"
     Test2  = suppressWarnings(wilcox.test(A ~ B))
     Z[k,10] = round(Test2$p.value, pDigits)
    }
   } # End Numeric and Binary
   
   ### START Ordinal and Binary
   
   if( Y[i,3] == "Binary" & Y[j,3] == "Ordinal" ){
     Z[k,5] = "Glass rank biserial"
     Test = suppressWarnings(wilcoxonRG(x = B, g = A, ci=ci, R=R, reportIncomplete=reportIncomplete, conf=conf))
     if(!ci){
      Z[k,6] = round(-1 * Test, digits)
     }
     if(ci){
      Z[k,6] = round(-1 * Test$rg, digits)
      Z[k,7] = round(-1 * Test$lower.ci, digits)
      Z[k,8] = round(-1 * Test$upper.ci, digits)
     }
     Z[k,9] = "wilcox.test"
     Test2  = suppressWarnings(wilcox.test(as.numeric(B) ~ A))
     Z[k,10] = round(Test2$p.value, pDigits)
    }
   
   if( Y[i,3] == "Ordinal" & Y[j,3] == "Binary" ){
    Z[k,5] = "Glass rank biserial"
    Test = suppressWarnings(wilcoxonRG(x = B, g = A, ci=ci, R=R, reportIncomplete=reportIncomplete, conf=conf))
    if(!ci){
     Z[k,6] = round(-1 * Test, digits)
    }
    if(ci){
     Z[k,6] = round(-1 * Test$rg, digits)
     Z[k,7] = round(-1 * Test$lower.ci, digits)
     Z[k,8] = round(-1 * Test$upper.ci, digits)
    }
    Z[k,9] = "wilcox.test"
    Test2  = suppressWarnings(wilcox.test(as.numeric(A) ~ B))
    Z[k,10] = round(Test2$p.value, pDigits)
   } 
   # End Ordinal and Binary
   
  ### START Numeric and Nominal
  
  if( Y[i,3] == "Nominal" & Y[j,3] == "Numeric" ){
    if(methodNumNom == "eta"){
     Z[k,5] = "Eta"
     Z[k,9] = "Anova"
     }
    if(methodNumNom == "epsilon"){
     Z[k,5] = "Epsilon"
     Z[k,9] = "Anova on ranks"
     B = rank(B)
    }
    model  = suppressWarnings(lm(B ~ A))
    Test   = summary(model)$r.squared
    Z[k,6] = round(sqrt(Test), digits)
    Params = length(levels(A)) - 1
    OBS    = length(A)
    SE     = sqrt((4 * Test * (1 - Test)^2 * (OBS - Params - 1)^2)/((OBS^2 - 1) * (OBS + 3)))
    Zvalue = -qnorm((1-conf)/2)
    lower.ci = Test - (Zvalue * SE)
      if(lower.ci < 0){lower.ci=0}
    upper.ci = Test + (Zvalue * SE)
    if(upper.ci > 1){upper.ci=1}
    Test2 = summary(model)$fstatistic
    Z[k,7] = round(sqrt(lower.ci), digits)
    Z[k,8] = round(sqrt(upper.ci), digits)
    Z[k,10] = round(pf(Test2[1],Test2[2],Test2[3],lower.tail=F), pDigits)
  }
   
   if( Y[i,3] == "Numeric" & Y[j,3] == "Nominal" ){
    if(methodNumNom == "eta"){
     Z[k,5] = "Eta"
     Z[k,9] = "Anova"
    }
    if(methodNumNom == "epsilon"){
     Z[k,5] = "Epsilon"
     Z[k,9] = "Anova on ranks"
     A = rank(A)
    }
    model  = suppressWarnings(lm(A ~ B))
    Test   = summary(model)$r.squared
    Z[k,6] = round(sqrt(Test), digits)
    Params = length(levels(A)) - 1
    OBS    = length(A)
    SE     = sqrt((4 * Test * (1 - Test)^2 * (OBS - Params - 1)^2)/((OBS^2 - 1) * (OBS + 3)))
    Zvalue = -qnorm((1-conf)/2)
    lower.ci = Test - (Zvalue * SE)
    if(lower.ci < 0){lower.ci=0}
    upper.ci = Test + (Zvalue * SE)
    if(upper.ci > 1){upper.ci=1}
    Test2 = summary(model)$fstatistic
    Z[k,7] = round(sqrt(lower.ci), digits)
    Z[k,8] = round(sqrt(upper.ci), digits)
    Z[k,10] = round(pf(Test2[1],Test2[2],Test2[3],lower.tail=F), pDigits)
   }

 } # End j for loop
 } # End i for loop

 if(ci){if(progress){cat("\n")}}
   
 for(h in 1:n){
  Z[h,"Signif"] = "n.s."
  
  if(is.na(Z[h,"p.value"])){Z[h,"Signif"] = "NA"}
  if(!is.na(Z[h,"p.value"])){   
     if( Z[h,"p.value"] <= 0.05 )   { Z[h,"Signif"] = "*"}
     if( Z[h,"p.value"] <= 0.01 )   { Z[h,"Signif"] = "**"}
     if( Z[h,"p.value"] <= 0.001 )  { Z[h,"Signif"] = "***"}
     if( Z[h,"p.value"] <= 0.0001 ) { Z[h,"Signif"] = "****"}
    }
 }
  
 return(Z)       
 }