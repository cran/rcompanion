#' @title Compact letter display for lists of comparisons
#'
#' @description Produces a compact letter display (cld) from pairwise 
#'              comparisons that were summarized in a table of comparisons
#'
#' @param formula A formula indicating the variable holding p-values and
#'                the variable holding the comparisons. 
#'                e.g. P.adj ~ Comparison.
#' @param data   The data frame to use.
#' @param comparison A vector of text describing comparisons, 
#'                   with each element in a form similar to
#'                   "Treat.A - Treat.B = 0".  Spaces and "=" and "0"
#'                   are removed by default
#' @param p.value    A vector of p-values corresponding to the comparisons
#'                   in the \code{comparison} argument
#' @param threshold  The alpha value.  That is, the p-value below which the
#'                   comparison will be considered significant
#' @param print.comp If \code{TRUE}, prints out a data frame of the
#'                   modified text of the comparisons.  Useful for debugging
#' @param remove.space  If \code{TRUE}, removes spaces from the text of the
#'                      comparisons
#' @param remove.equal  If \code{TRUE}, removes "=" from the text of the
#'                      comparisons
#' @param remove.zero   If \code{TRUE}, removes "0" from the text of the
#'                      comparisons
#' @param swap.colon    If \code{TRUE}, swaps ":" with "-" in the text of the
#'                      comparisons
#' @param swap.vs       If \code{TRUE}, swaps "vs" with "-" in the text of the
#'                      comparisons                      
#' @param ...           Additional arguments passed to
#'                      \code{multcompLetters}              
#'             
#' @details  The input should include either \code{formula} and \code{data};
#'           or \code{comparison} and \code{p.value}.
#'           
#'           This function relies upon the \code{multcompLetters}
#'           function in the \code{multcompView} package. The text for the
#'           comparisons
#'           passed to \code{multcompLetters} should be in the form
#'           "Treat.A-Treat.B".  Currently by default \code{cldList} removes
#'           spaces, equal signs, and zeros, by default, 
#'           and so can use 
#'           text in the form e.g.
#'           "Treat.A - Treat.B = 0".
#'           It also changes ":" to "-", and so can use
#'           text in the form e.g.
#'           "Treat.A : Treat.B".
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_08.html}
#' @concept compact letter display cld post-hoc
#' @return A data frame of group names, group separation letters,
#'         and monospaced separtions letters
#'
#' @note  The parsing of the formula is simplistic. 
#'          The first variable on the
#'          left side is used as the measurement variable.  
#'          The first variable on the
#'          right side is used for the grouping variable.  
#'      
#'       It is often helpful to reorder the factor levels in the data
#'       set so that the group with the largest e.g.
#'       mean or median is first, and so on.
#'      
#' @examples
#' data(BrendonSmall)
#' 
#' model = aov(Calories ~ Instructor, data=BrendonSmall)
#' 
#' TUK = TukeyHSD(model, "Instructor", ordered = TRUE)
#' 
#' ### Convert the TukeyHSD output to a standard data frame
#' 
#' TUK = as.data.frame(TUK$Instructor)
#' names(TUK) = gsub(" ", ".", names(TUK))
#' 
#' HSD = data.frame(Comparison=row.names(TUK), 
#'                  diff=TUK$diff, lwr=TUK$lwr, lwr=TUK$lwr, p.adj=TUK$p.adj)
#' 
#' HSD
#' 
#' cldList(p.adj ~ Comparison, data = HSD,
#'         threshold = 0.05,
#'         remove.space=FALSE)
#' 
#' @importFrom multcompView multcompLetters
#' 
#' @export
 
cldList = function(formula       = NULL,
                   data          = NULL, 
                   comparison    = NULL, 
                   p.value       = NULL, 
                   threshold     = 0.05,
                   print.comp    = FALSE,
                   remove.space  = TRUE,
                   remove.equal  = TRUE,
                   remove.zero   = TRUE,
                   swap.colon    = TRUE,
                   swap.vs       = FALSE,
                   ...)
{
  if(!is.null(formula)){
    p.value     = eval(parse(text=paste0("data","$",all.vars(formula[[2]])[1])))
    comparison  = eval(parse(text=paste0("data","$",all.vars(formula[[3]])[1])))
  }
  
FLAG = 0
  
Comparison = (as.numeric(p.value) <= threshold)

if (sum(Comparison) == 0){FLAG =1}

if(remove.space == TRUE) {comparison = gsub(" ",  "",  comparison)}
if(remove.equal == TRUE) {comparison = gsub("=",  "",  comparison)}
if(remove.zero  == TRUE) {comparison = gsub("0",  "",  comparison)}
if(swap.colon   == TRUE) {comparison = gsub(":",  "-", comparison)}
if(swap.vs      == TRUE) {comparison = gsub("vs", "-", comparison)}

names(Comparison) = comparison

if(print.comp == TRUE) 
  {Y = data.frame(Comparisons = names(Comparison), 
                  p.value = p.value, Value=Comparison,
                  Threshold=threshold)
   cat("\n", "\n")
   print(Y)
   cat("\n", "\n")}

MCL = multcompLetters(Comparison, ...)

Group      = names(MCL$Letters)
Letter     = as.character(MCL$Letters)
if(FLAG==0){MonoLetter = as.character(MCL$monospacedLetters)}
if(FLAG==1){MonoLetter = Letter}

Z = data.frame(Group, Letter, MonoLetter)

return(Z)
}