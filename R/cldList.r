#' @title Compact letter display for lists of comparisons
#'
#' @description Produces a compact letter display (cld) from pairwise 
#'              comparisons that were summarized in a table of comparisons
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
#' @param ...           Additional arguments passed to
#'                      \code{multcompLetters}              
#'             
#' @details  This function relies upon the \code{multcompLetters}
#'           function in the \code{multcompView} package. The text for the
#'           comparisons
#'           passed to \code{multcompLetters} should be in the form
#'           "Treat.A-Treat.B".  Currently \code{cldList} removes
#'           spaces, equal signs, and zeros by default, and so can use 
#'           text in the form e.g.
#'           "Treat.A - Treat.B = 0".
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' @references \url{http://rcompanion.org/handbook/F_08.html}
#' @concept compact letter display cld post-hoc
#' @return A data frame of group names, group separation letters,
#'         and monospaced separtions letters
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
#' cldList(comparison       = DT$Comparison,
#'         p.value          = DT$P.adj,
#'         threshold        = 0.05)
#' 
#' @importFrom multcompView multcompLetters
#' 
#' @export
 
cldList = function(comparison, 
                   p.value, 
                   threshold     = 0.05,
                   print.comp    = FALSE,
                   remove.space  = TRUE,
                   remove.equal  = TRUE,
                   remove.zero   = TRUE,
                   ...)
{
Comparison = p.value < threshold

if(remove.space == TRUE) {comparison = gsub(" ", "", comparison)}
if(remove.equal == TRUE) {comparison = gsub("=", "", comparison)}
if(remove.zero  == TRUE) {comparison = gsub("0", "", comparison)}

names(Comparison) = comparison

if(print.comp == TRUE) 
  {Y = data.frame(Comparisons = names(Comparison))
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