#' @title Functions to Support Extension Education Program Evaluation
#' 
#' @description Functions and datasets to support Summary and Analysis of
#' Extension Program Evaluation in R and An R
#' Companion for the Handbook of Biological Statistics.
#'
#' @section  Useful functions:
#'
#' There are several functions that provide summary statistics for
#' grouped data. These function titles tend to start with \code{"groupwise"}.
#' They provide means, medians, geometric means, and Huber M-estimators
#' for groups, along with confidence intervals by traditional
#' methods and bootstrap.
#' 
#' Functions to produce effect size statistics, 
#' some with bootstrapped confidence intervals, 
#' include those for Cramer's V, 
#' Cohen's g and odds ratio for paired tables, Cohen's h, 
#' Cohen's w, Vargha and Delaney's A, Cliff's delta,
#' r for one-sample, two-sample, and paired Wilcoxon and Mann-Whitney tests, 
#' epsilon-squared, and Freeman's theta. 
#'
#' There are also functions that are useful for comparing models.
#' \code{\link{compareLM}},  \code{\link{compareGLM}}, and 
#' \code{\link{pairwiseModelAnova}}.
#' These use goodness-of-fit measures like AIC, BIC, and BICc, or likelihood
#' ratio tests. The \code{\link{accuracy}} function reports statistics for
#' models including minimum maximum accuracy, MAPE, RMSE, 
#' Efron's pseudo r-squared, and coefficient of variation.
#'
#' Functions for nominal data include post-hoc tests for 
#' Cochran-Mantel-Haenszel test (\code{\link{groupwiseCMH}}),
#' for McNemar-Bowker test (\code{\link{pairwiseMcnemar}}),
#' and for tests of association like Chi-square, Fisher exact, and G-test
#' (\code{\link{pairwiseNominalIndependence}}).
#' 
#' There are a few useful plotting functions, including 
#' \code{\link{plotNormalHistogram}} that plots a histogram of values and 
#' overlays
#' a normal curve, and \code{\link{plotPredy}} which plots of line for predicted
#' values for a bivariate model.  Other plotting functions include producing
#' density plots.
#'
#' The function \code{\link{nagelkerke}}
#' provides pseudo R-squared values for a variety of model types, as well as
#' a likelihood ratio test for the model as a whole.
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
#' @docType package
#' 
#' @name rcompanion
#' 
NULL