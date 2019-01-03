#' @section  Useful functions:
#' The function \code{\link{nagelkerke}}
#' provides pseudo R-squared values for a variety of model types, as well as
#' a likelihood ratio test for the model as a whole.
#'
#' There are several functions that provide summary statistics for
#' grouped data. These function titles tend to start with \code{"groupwise"}.
#' They provide means, medians, geometric means, and Huber M-estimators
#' for groups, along with confidence intervals by traditional
#' methods and bootstrap.
#'
#' Function titles starting with \code{"pairwise"} conduct pairwise
#' tests among groups as a post-hoc analysis for omnibus tests.
#' These tests are Mood's median test, sign test (for omnibus Friedman test),
#' and permutation test.
#' The output can be parsed into a compact letter display.
#'
#' There are also functions that are useful for comparing models.
#' \code{\link{compareLM}},  \code{\link{compareGLM}}, and 
#' \code{\link{pairwiseModelAnova}}.
#' These use goodness-of-fit measures like AIC, BIC, and BICc, or likelihood
#' ratio tests.
#'
#' There are a few useful plotting functions, including 
#' \code{\link{plotNormalHistogram}} that plots a histogram of values and 
#' overlays
#' a normal curve, and \code{\link{plotPredy}} which plots of line for predicted
#' values for a bivariate model.  Other plotting functions include producing
#' density plots.
#'
#' Functions for nominal data include post-hoc tests for 
#' Cochran-Mantel-Haenszel test (\code{\link{groupwiseCMH}}),
#' for McNemar-Bowker test (\code{\link{pairwiseMcnemar}}),
#' and for tests of association like Chi-square, Fisher exact, and G-test
#' (\code{\link{pairwiseNominalIndependence}}).
#'
#' A function close to my heart is (\code{\link{cateNelson}}), which performs
#' Cate-Nelson analysis for bivariate data.
#'
#' @section Vignettes and examples:
#' The functions in this package are used in
#' "Extension Education Program Evaluation in R" which is available at
#' \url{http://rcompanion.org/handbook/}
#' and "An R Companion for the Handbook of Biological Statistics" 
#' which is available at \url{http://rcompanion.org/rcompanion/}.
#' 
#' The documentation for each function includes an example as well.
#'
#' @section  Version notes:
#' Version 2.0 is not entirely back-compatable
#' as several functions have been removed.
#' These include some of the pairwise methods that can be replaced with
#' better methods.  Also, some functions have been removed or modified
#' in order to import fewer packages.
#' 
#' Removed packages are indicated with 'Defunct' in their titles.
#'
#'@keywords internal
"_PACKAGE"