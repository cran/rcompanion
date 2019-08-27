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
#'           "Treat.A-Treat.B".  Currently \code{cldList} removes
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
#' @examples
#' data(PoohPiglet)
#' PoohPiglet$Speaker = factor(PoohPiglet$Speaker,
#'                             levels=c("Pooh", "Tigger", "Piglet"))
#' library(FSA)
#' DT = dunnTest(Likert ~ Speaker,
#'               data=PoohPiglet,
#'               method="bh")
#' DT = DT$res
#' DT
#' cldList(P.adj ~ Comparison,
#'         data      = DT,
#'         threshold = 0.05)
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
Comparison = (as.numeric(p.value) <= threshold)

if (sum(Comparison) == 0){stop("No significant differences.", call.=FALSE)}

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
MonoLetter = as.character(MCL$monospacedLetters)

Z = data.frame(Group, Letter, MonoLetter)

return(Z)
}