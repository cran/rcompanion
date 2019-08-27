#' @title Wilcoxon z statistic
#' 
#' @description Calculates the z statistic for a Wilcoxon
#'              two-sample, paired, or one-sample test.
#' 
#' @param x A vector of observations.
#' @param y For the two-sample and paired cases,
#'          a second vector of observations.
#' @param mu For the one-sample case,
#'           the value to compare \code{x} to, as in \code{wilcox.test}
#' @param paired  As used in \code{wilcox.test}.
#' @param exact   As used in \code{wilcox.test}, 
#'                default here is \code{FALSE}.
#' @param correct As used in \code{wilcox.test}, 
#'                default here is \code{FALSE}.              
#' @param digits  The number of significant digits in the output.
#'
#' @details  This function uses code from \code{wilcox.test},
#'           and reports the \code{z} statistic,
#'           which is calculated by the original function
#'           but isn't returned.
#'           
#'           The returned value will be NA if the function attempts an
#'           exact test.
#'           
#'           For the paired case, the observations in \code{x} and
#'           and \code{y} should be ordered such that the
#'           first observation in \code{x} is paired with the first observation
#'           in \code{y}, and so on.
#'           
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}, 
#'         R Core Team
#' @return A single statistic, \code{z}.
#'         
#' @examples
#' data(Pooh)
#' wilcoxonZ(x = Pooh$Likert[Pooh$Time==1], y = Pooh$Likert[Pooh$Time==2],
#'           paired=TRUE, exact=FALSE, correct=FALSE)
#' 
#' @importFrom stats setNames
#' 
#' @export

wilcoxonZ <-
function(x, y=NULL,
         mu=0, paired=FALSE, exact=FALSE, correct=FALSE,
         digits=3)
 {
    if(!missing(mu) && ((length(mu) > 1L) || !is.finite(mu)))
        stop("'mu' must be a single number")
    if(!is.numeric(x)) stop("'x' must be numeric")
    if(!is.null(y)) {
        if(!is.numeric(y)) stop("'y' must be numeric")
        if(paired) {
            if(length(x) != length(y))
                stop("'x' and 'y' must have the same length")
            OK <- complete.cases(x, y)
            x <- x[OK] - y[OK]
            y <- NULL
        }
        else {
            x <- x[is.finite(x)]
            y <- y[is.finite(y)]
        }
    } else {
        if(paired)
            stop("'y' is missing for paired test")
        x <- x[is.finite(x)]
    }

    if(length(x) < 1L)
        stop("not enough (finite) 'x' observations")
    CORRECTION <- 0
    if(is.null(y)) {
        x <- x - mu
        ZEROES <- any(x == 0)
        if(ZEROES)
            x <- x[x != 0]
        n <- as.double(length(x))
        if(is.null(exact)){exact <- (n < 50)}
        r <- rank(abs(x))
        STATISTIC <- setNames(sum(r[x > 0]), "V")
        TIES <- length(r) != length(unique(r))

        if(exact && !TIES && !ZEROES) {
          z = NA
        } else { ## not exact, maybe ties or zeroes
            NTIES <- table(r)
            z <- STATISTIC - n * (n + 1)/4
            SIGMA <- sqrt(n * (n + 1) * (2 * n + 1) / 24
                          - sum(NTIES^3 - NTIES) / 48)
            if(correct) {
                CORRECTION <-sign(z) * 0.5
            }
            z <- (z - CORRECTION) / SIGMA
        }
    }
        else { ##-------------------------- 2-sample case ---------------------
        if(length(y) < 1L)
            stop("not enough 'y' observations")
        r <- rank(c(x - mu, y))
        n.x <- as.double(length(x))
        n.y <- as.double(length(y))
        if(is.null(exact)){exact <- (n.x < 50) && (n.y < 50)}
        STATISTIC <- c("W" = sum(r[seq_along(x)]) - n.x * (n.x + 1) / 2)
        TIES <- (length(r) != length(unique(r)))
        if(exact && !TIES) {
          z = NA
            }
        else {
            NTIES <- table(r)
            z <- STATISTIC - n.x * n.y / 2
            SIGMA <- sqrt((n.x * n.y / 12) *
                          ((n.x + n.y + 1)
                           - sum(NTIES^3 - NTIES)
                           / ((n.x + n.y) * (n.x + n.y - 1))))
            if(correct) {CORRECTION <- sign(z) * 0.5}
	    z <- (z - CORRECTION) / SIGMA
      
            if(exact && TIES) {
                warning("cannot compute exact p-value with ties")
            }
        }
        }
	    z = signif(z, digits=digits)
	    names(z) = "z"
	    return(z)
}