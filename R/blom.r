#' @title Normal scores transformation
#'
#' @description Normal scores transformation
#'              (Inverse normal transformation)
#'              by Elfving, Blom, van der Waerden, Tukey, 
#'              and rankit methods,
#'              as well as z score transformation
#'              (standardization)
#'              and scaling to a range (normalization).
#' 
#' @param x A vector of numeric values.
#' @param method Any one \code{"general"} (the default), 
#'               \code{"blom"}, \code{vdw},
#'               \code{"tukey"}, \code{"elfving"},
#'               \code{"rankit"},
#'               \code{zscore}, or \code{scale}.
#' @param alpha A value used in the \code{"general"} method.
#'          If alpha=pi/8 (the default), the \code{"general"} method reduces  
#'          to the \code{"elfving"} method.
#'          If alpha=3/8, the \code{"general"} method reduces 
#'          to the \code{"blom"} method.
#'          If alpha=1/2, the \code{"general"} method reduces to the 
#'          \code{"rankit"} method.
#'          If alpha=1/3, the \code{"general"} method reduces to the 
#'          \code{"tukey"} method.
#'          If alpha=0, the \code{"general"} method reduces to the 
#'          \code{"vdw"} method.
#' @param complete  If \code{TRUE}, \code{NA} values are removed
#'                  before transformation. The default is \code{FALSE}.
#' @param na.last Passed to \code{rank} in the normal scores methods. 
#'                See the documentation for the
#'                \code{rank} function.
#'                The default is \code{"keep"}.
#' @param na.rm  Used in the \code{"zscore"} and \code{"scale"} methods.
#'               Passed to \code{mean}, \code{min}, and \code{max}
#'               functions in those methods.
#'               The default is \code{TRUE}.
#' @param adjustN  If \code{TRUE}, the default, the normal scores methods 
#'                 use only non-\code{NA} values to determine the sample size,
#'                 \code{N}.  This seems to work well under default conditions
#'                 where \code{NA} values are retained, even if there are
#'                 a high percentage of \code{NA} values.
#' @param min  For the \code{"scale"} method, the minimum value of the
#'             transformed values.
#' @param max  For the \code{"scale"} method, the maximum value of the
#'             transformed values.
#' @param ...  additional arguments passed to \code{rank}.
#' 
#' @details By default, \code{NA} values are retained in the output.
#'          This behavior can be changed with the \code{na.rm} argument
#'          for \code{"zscore"} and \code{"scale"} methods, or
#'          with \code{na.last} for the normal scores methods.
#'          Or \code{NA} values can be removed from the input with
#'          \code{complete=TRUE}.
#'          
#'          For normal scores methods, if there are \code{NA} values
#'          or tied values,
#'          it is helpful to look up 
#'          the documentation for \code{rank}.
#'          
#'          In general, for normal scores methods, either of the arguments
#'          \code{method} or \code{alpha} can be used.
#'          With the current algorithms, there is no need to use both.
#'          
#'          Normal scores transformation will return a
#'          normal distribution with a mean of 0 and
#'          a standard deviation of 1.
#'          
#'          The \code{"scale"} method coverts values to the range specified
#'          in \code{max} and \code{min} without transforming the distribution
#'          of values. By default, the \code{"scale"} method converts values
#'          to a 1 to 10 range.
#'          Using the \code{"scale"} method with
#'          \code{min = 0} and \code{max = 1} is
#'          sometimes called "normalization".
#'          
#'          The \code{"zscore"} method converts values by the usual method
#'          for z scores: \code{(x - mean(x)) / sd(x)}. The transformed
#'          values with have a mean of 0 and a standard deviation of
#'          1 but won't be coerced into a normal distribution.
#'          Sometimes this method is called "standardization".
#'          
#' @author Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
#' 
#' @references Conover, 1995, Practical Nonparametric Statistics, 3rd.
#' 
#'             Solomon & Sawilowsky, 2009, 
#'             Impact of rank-based normalizing transformations 
#'             on the accuracy of test scores.
#'             
#'             Beasley and Erickson, 2009, Rank-based inverse normal 
#'             transformations are increasingly used, but are they merited?
#'
#' @note It's possible that Gustav Elfving didn't recommend the 
#'       formula used in this function for the \code{Elfving} method.  
#'       I would like thank Terence Cooke
#'       at the University of Exeter for their
#'       diligence at trying to track down a reference for this formula. 
#'
#' @concept normal scores
#' @concept z score
#' @concept standardization
#' @concept normalization
#' @concept van der Waerden
#' @concept Blom
#' @concept Elfving
#' @concept rankit
#' 
#' @return A vector of numeric values.
#'           
#' @examples
#' set.seed(12345)
#' A = rlnorm(100)
#' \dontrun{hist(A)}
#' ### Convert data to normal scores by Elfving method
#' B = blom(A)
#' \dontrun{hist(B)}
#' ### Convert data to z scores 
#' C = blom(A, method="zscore")
#' \dontrun{hist(C)}
#' ### Convert data to a scale of 1 to 10 
#' D = blom(A, method="scale")
#' \dontrun{hist(D)}
#' 
#' ### Data from Sokal and Rohlf, 1995, 
#' ### Biometry: The Principles and Practice of Statistics
#' ### in Biological Research
#' Value = c(709,679,699,657,594,677,592,538,476,508,505,539)
#' Sex   = c(rep("Male",3), rep("Female",3), rep("Male",3), rep("Female",3))
#' Fat   = c(rep("Fresh", 6), rep("Rancid", 6))
#' ValueBlom = blom(Value)
#' Sokal = data.frame(ValueBlom, Sex, Fat)
#' model = lm(ValueBlom ~ Sex * Fat, data=Sokal)
#' anova(model)
#' \dontrun{
#' hist(residuals(model))
#' plot(predict(model), residuals(model))
#' }
#' 
#' @importFrom stats qnorm
#' 
#' @export

blom = function(x, method="general", alpha=pi/8, 
                complete=FALSE, na.last="keep", na.rm=TRUE,
                adjustN=TRUE,
                min=1, max=10, ...){
  if(complete){x=x[complete.cases(x)]}
  Ranks = rank(x, na.last=na.last, ...)
  if(adjustN==FALSE){N = length(x)}
  if(adjustN==TRUE) {N = sum(complete.cases(x))}
  if(method=="blom")   {Score = qnorm((Ranks-0.375)/(N+0.25))}
  if(method=="vdw")    {Score = qnorm((Ranks)/(N+1))}
  if(method=="tukey")  {Score = qnorm((Ranks-1/3)/(N+1/3))}
  if(method=="rankit") {Score = qnorm((Ranks-1/2)/(N))}
  if(method=="elfving"){Score = qnorm((Ranks-pi/8)/(N-pi/4+1))}
  if(method=="general"){Score = qnorm((Ranks-alpha)/(N-2*alpha+1))}
  if(method=="zscore") {Score = (x-mean(x, na.rm=na.rm))/sd(x, na.rm=na.rm)}
  if(method=="scale")  {
    Score = (((x - min(x, na.rm=na.rm)) * 
            (max - min)) / 
            (max(x, na.rm=na.rm) - min(x, na.rm=na.rm)) + 
             min)
  }
  return(Score)
}