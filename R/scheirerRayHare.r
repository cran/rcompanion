#' @title Scheirer Ray Hare test
#'
#' @description Conducts Scheirer Ray Hare test.
#'              
#' @param formula A formula indicating the response variable and
#'                two independent variables. e.g. y ~ x1 + x2.
#' @param data   The data frame to use. 
#' @param y If no formula is given, the response variable.
#' @param x1 If no formula is given, the first independent variable.
#' @param x2 If no formula is given, the second independent variable.
#' @param tie.correct If \code{"TRUE"}, applies a correction for ties in the 
#'                                      response variable.
#' @param ss If \code{"TRUE"}, includes the sums of squares in the output.
#' @param verbose If \code{"TRUE"}, outputs statistics used in the analysis 
#'                                  by direct print.
#'
#' @details The Scheirer Ray Hare test is a nonparametric test used for a 
#'          two-way factorial experiment.  It is described by Sokal and
#'          Rohlf (1995).
#'          It is sometimes recommended that the design should be balanced,
#'          and that there should be at least five observations for each
#'          cell in the interaction.
#'          One might consider using aligned ranks transformation anova
#'          instead of the Scheirer Ray Hare test.
#'          
#'          The input should include either \code{formula} and \code{data};
#'          or \code{y}, \code{x1}, and \code{x2}.
#'          
#'          The function removes cases with NA in any of the variables.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references Sokal, R.R. and F.J. Rohlf. 1995. Biometry. 3rd ed. W.H. Freeman, 
#'             New York.
#'             
#'             \url{http://rcompanion.org/handbook/F_14.html}
#'             
#' @concept nonparametric Scheirer Ray Hare two-way
#' @return A data frame of results similar to an anova table. Output from the
#'         \code{verbose} option is printed directly and not returned with
#'         the data frame.
#'         
#' @note    The parsing of the formula is simplistic. 
#'          The first variable on the
#'          left side is used as the measurement variable.  
#'          The first variable on the
#'          right side is used for the first independent variable.
#'          The second variable on the
#'          right side is used for the second independent variable.  
#'  
#'         
#' @examples
#' ### Example from Sokal and Rohlf, 1995.
#' Value = c(709,679,699,657,594,677,592,538,476,508,505,539)
#' Sex   = c(rep("Male",3), rep("Female",3), rep("Male",3), rep("Female",3))
#' Fat   = c(rep("Fresh", 6), rep("Rancid", 6))
#' Sokal = data.frame(Value, Sex, Fat)
#' 
#' scheirerRayHare(Value ~ Sex + Fat, data=Sokal)
#'                                                              
#' @importFrom utils tail
#' @importFrom stats lm pchisq anova complete.cases
#' 
#' @export

scheirerRayHare = 
  function(formula=NULL, data=NULL,
           y=NULL, x1=NULL, x2=NULL, 
           tie.correct=TRUE, ss=TRUE, verbose=TRUE)
  {
  if(is.null(formula)){
    yname  = tail(as.character(substitute(y)), n=1)
    x1name = tail(as.character(substitute(x1)), n=1)
    x2name = tail(as.character(substitute(x2)), n=1)
  }
  if(!is.null(formula)){
    yname  = all.vars(formula[[2]])[1]
    x1name = all.vars(formula[[3]])[1]
    x2name = all.vars(formula[[3]])[2]
    y   = eval(parse(text=paste0("data","$",yname)))
    x1  = eval(parse(text=paste0("data","$",x1name)))
    x2  = eval(parse(text=paste0("data","$",x2name)))
    }
    
  Complete = complete.cases(y, x1, x2)
  y  = y[Complete]
  x1 = x1[Complete]
  x2 = x2[Complete]
  if(!is.factor(x1)){x1=factor(x1)}
  if(!is.factor(x2)){x2=factor(x2)} 
 
  Ranks = rank(y)
  Ties  = table(Ranks)
  Model = lm(Ranks ~ x1 + x2 + x1:x2)
  Anva  = anova(Model)
  MS    = Anva[1:4,1:3]
  n     = length(Ranks)
  D     = 1
  if(tie.correct){D = (1 - sum(Ties^3 - Ties)/(n^3 - n))}
  MStotalSokal   = n*(n+1)/12
  SS      = MS[1:3,2]
  # SStotal = sum(MS[,2])
  # DFtotal = sum(MS[,1])
  # MStotal = SStotal / DFtotal
  H = SS / MStotalSokal
  Hadj = H / D
  MS[1:3,4] =  Hadj
  MS[1:3,5] = (1-pchisq(MS[1:3,4],MS[1:3,1]))
  
    if(verbose){cat("\n")
              cat("DV: ", yname, "\n")
              cat("Observations: ", n, "\n")
              cat("D: ", D, "\n")
              cat("MS total: ", MStotalSokal, "\n")
              cat("\n")
              }
  
  colnames(MS)[4:5] = c("H","p.value")
  rownames(MS) = c(x1name, x2name, paste0(x1name,":",x2name), "Residuals")

 if(!ss){Z = MS[,c(1,4,5)]}
 if(ss){Z  = MS[,c(1,2,4,5)]}
 return(Z)
}