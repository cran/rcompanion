% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scheirerRayHare.r
\name{scheirerRayHare}
\alias{scheirerRayHare}
\title{Scheirer Ray Hare test}
\usage{
scheirerRayHare(
  formula = NULL,
  data = NULL,
  y = NULL,
  x1 = NULL,
  x2 = NULL,
  type = 2,
  tie.correct = TRUE,
  ss = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{formula}{A formula indicating the response variable and
two independent variables. e.g. y ~ x1 + x2.}

\item{data}{The data frame to use.}

\item{y}{If no formula is given, the response variable.}

\item{x1}{If no formula is given, the first independent variable.}

\item{x2}{If no formula is given, the second independent variable.}

\item{type}{The type of sum of squares to be used.
Acceptable options are \code{1}, \code{2},
\code{"I"}, or \code{"II"}.}

\item{tie.correct}{If \code{"TRUE"}, applies a correction for ties in the 
response variable.}

\item{ss}{If \code{"TRUE"}, includes the sums of squares in the output.}

\item{verbose}{If \code{"TRUE"}, outputs statistics used in the analysis 
by direct print.}
}
\value{
A data frame of results similar to an anova table. Output from the
        \code{verbose} option is printed directly and not returned with
        the data frame.
}
\description{
Conducts Scheirer Ray Hare test.
}
\details{
The Scheirer Ray Hare test is a nonparametric test used for a 
         two-way factorial experiment.  It is described by Sokal and
         Rohlf (1995).
         
         It is sometimes recommended that the design should be balanced,
         and that there should be at least five observations for each
         cell in the interaction.
         
         One might consider using aligned ranks transformation anova
         instead of the Scheirer Ray Hare test.
         
         Note that for unbalanced designs, 
         by default, a type-II sum-of-squares
         approach is used.
         
         The input should include either \code{formula} and \code{data};
         or \code{y}, \code{x1}, and \code{x2}.
         
         The function removes cases with NA in any of the variables.
}
\note{
The parsing of the formula is simplistic. 
         The first variable on the
         left side is used as the measurement variable.  
         The first variable on the
         right side is used for the first independent variable.
         The second variable on the
         right side is used for the second independent variable.
}
\section{Acknowledgments}{

 
         Thanks to Guillaume Loignon for the suggestion to
         include type-II sum-of-squares.
}

\examples{
### Example from Sokal and Rohlf, 1995.
Value = c(709,679,699,657,594,677,592,538,476,508,505,539)
Sex   = c(rep("Male",3), rep("Female",3), rep("Male",3), rep("Female",3))
Fat   = c(rep("Fresh", 6), rep("Rancid", 6))
Sokal = data.frame(Value, Sex, Fat)

scheirerRayHare(Value ~ Sex + Fat, data=Sokal)
                                                             
}
\references{
Sokal, R.R. and F.J. Rohlf. 1995. Biometry. 3rd ed. W.H. Freeman, 
            New York.
            
            \url{https://rcompanion.org/handbook/F_14.html}
}
\author{
Salvatore Mangiafico, \email{mangiafico@njaes.rutgers.edu}
}
\concept{Scheirer-Ray-Hare}
