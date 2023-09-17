#' @title Convert a lower triangle matrix to a full matrix
#'
#' @description Converts a lower triangle matrix to a full matrix. 
#' 
#' @param PT A lower triangle matrix.
#' 
#' @details  This function is useful to convert a lower triangle matrix
#'           of p-values from a pairwise test to a full matrix.
#'           A full matrix can be passed to \code{multcompLetters}
#'           in the \code{multcompView} package to produce a compact
#'           letter display.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references \url{https://rcompanion.org/handbook/F_08.html}
#' 
#' @seealso \code{\link{cldList}}
#' 
#' @concept multiple comparisons
#' @concept compact letter display
#' @concept cld
#' @concept post-hoc
#' 
#' @return A full matrix.
#'         
#' @examples
#' ### Example with pairwise.wilcox.test
#' data(BrendonSmall)
#' BrendonSmall$Instructor = factor(BrendonSmall$Instructor,
#'                           levels = c('Brendon Small', 'Jason Penopolis',
#'                                      'Paula Small', 'Melissa Robbins', 
#'                                      'Coach McGuirk'))
#' P   = pairwise.wilcox.test(x = BrendonSmall$Score, g = BrendonSmall$Instructor)
#' PT  = P$p.value
#' PT
#' PT1 = fullPTable(PT)
#' PT1
#' library(multcompView)
#' multcompLetters(PT1)
#' 
#' @importFrom multcompView multcompLetters
#' 
#' @export

fullPTable = 
  function(PT)
  {
PTa <- rep(c(NA_real_),length(PT[,1]))         # Add a row
PT1 <- rbind(PTa, PT)
rownames(PT1)[1] <- colnames(PT1)[1]
PTb <- rep(c(NA_real_),length(PT1[,1]))        # Add a column
PT1 <- cbind(PT1, PTb)
n <- length(PT1[,1])
colnames(PT1)[n] <- rownames(PT1)[n]
PT2 = t(PT1)                                   # Create new transposed
PT2[lower.tri(PT2)] = PT1[lower.tri(PT1)]
diag(PT2) = signif(1.00, digits = 4)           # Set the diagonal
PT2
  }