#' rcompanion: Functions and data sets for rcompanion.org
#'
#'This package provides custom functions for working through examples 
#'and analyses from "Summary and Analysis of
#'Extension Education Program Evaluation in R" and "An R
#'Companion for the Handbook of Biological Statistics". 
#'             
#'@section  Useful functions:
#'The function \code{nagelkereke} 
#'provides pseudo R-squared values for variety of model types, as well as
#'a likelihood ratio test for the model as a whole.  An addtional function,
#'\code{nagelkerkeHermite}, is provided for models fit with the
#'\code{hermite} package.
#'
#'There are several functions which provide summary statistics for
#'grouped data. These function titles tend to start with \code{"groupwise"}.
#'They provide means, medians, geometric means, and Huber M-estimators
#'for groups, along with confidence intervals by traditional
#'methods and bootstrap.
#'
#'Function titles starting with \code{"pairwise"} conduct pairwise
#'tests among groups as a post-hoc analysis for omnibus tests.
#'At the time of writing, 
#'these tests are Mood's median test, sign test (for Friedman test),
#'permutation test, robust anova, and ordinal regression.
#'The output is a table of comparisons and p-values,
#'or a matrix of p-values that can be parsed into
#'a compact letter display.
#'
#'There are also functions that are useful for comparing models.
#'\code{compareLM},  \code{compareGLM}, and \code{pairwiseModelAnova}.
#'These use goodness-of-fit measures like AIC, BIC, and BICc, or likelihood
#'ratio tests.
#'
#'There are a few useful plotting functions, including 
#'\code{plotNormalHistogram} that plots a histogram of values and overlays
#'a normal curve, and \code{plotPredy} which plots of line for predicted
#'values for a bi-variate model.  Other plotting functions include producing
#'density plots.
#'
#'Functions for mominal data include post-hoc tests for 
#'Cochran-Mantel-Haenszel test (\code{groupwiseCMH}),
#'for McNemar-Bowker (\code{pairwiseMcneamr}),
#'and for tests of association like Chi-square, Fisher exact, and G-test
#'(\code{pairwiseNominalIndependence}).
#'
#'A function close to my heart is (\code{cateNelson}), which performs
#'Cate-Nelson analysis for bivariate data.
#'
#'@section Vignettes and examples:
#'The functions in this package are used in
#' "Extension Education Program Evaluation in R" which is available at
#' http://rcompanion.org/handbook/
#' and "An R Companion for the Handbook of Biological Statistics" 
#' which is available at http://rcompanion.org/rcompanion/.
#' 
#' The documentation for each function includes an example as well.
#'
#' @name rcompanion
#' @docType package
NULL
