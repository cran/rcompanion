## rcompanion v 2.4.36 (2024-05-27)

-   An option was added to the `pairwiseNominalIndependence()` function
    to allow for not using Yates' correction in 2 x 2 chi-square test.
    
## rcompanion v 2.4.35 (2024-02-15)

-   A fix was applied to the `transformTukey()` function
    to allow for samples sizes greater than 5000.
    
-   A fix was applied to the `cramerV()` function
    to not round the value for the total sample size.

## rcompanion v 2.4.34 (2023-09-15)

-   A fix was applied to the `cramerV()` function to avoid an error 
    in the case of certain sparse matrices.
    
-   An update was made in the example for `wilcoxonPairedR()` 
    to reflect the handling of
    paired samples in `wilcox.test()`. Specifically, this is in reference to 
    PR 14359, the correction for which doesn't allow paired arguments for 
    formula input for the `t.test()` and `wilcox.test()` functions.
    
- URLs in the documentation were updated to https where appropriate.

## rcompanion v 2.4.30 (2023-05-03)

-   Package documentation updated.

## rcompanion v 2.4.26 (2023-04-03)

-   MASS::rlm models were added to the accuracy() function.

-   The ephronRSquared() function now accepts model objects.

-   The accuracy(), ephronRSquared(), and nagelkerke() functions now accept glmmTMB model objects.

## rcompanion v 2.4.21 (2023-01-16)

-   Minor updates

## rcompanion v 2.4.18 (2022-08-06)

-   The `wilcoxonOneSampleR()` and `wilcoxonPairedR()` functions were modified to adjust the sample size for tied values.

-   Two functions, `mangiaficoD()` and `multiMangiaficoD()` were added. These functions calculate the difference in medians for two groups divided by pooled median absolute deviation. This is somewhat analogous to a nonparametric version of Cohen's d.

## rcompanion v 2.4.16 (2022-07-04)

-   The `pairwiseNominalIndependence()` function was modified to include chi-square, G, and df values in the output.

## rcompanion v 2.4.15 (2022-03-14)

-   The `scheirerRayHare()` function was modified to include an option to use type-II sum-of-squares.

-   The `cliffDelta()` function has been modified to better handle large sample sizes.

## rcompanion v 2.4.13 (2022-01-03)

-   The `cateNelsonFixedY()` function was added. This conducts a Cate-Nelson analysis on bivariate data where the critical Y value is fixed. This kind of analysis is useful in agronomy studies (for example, when a 95% yield is of interest), and for other bivariate data which don't conform well to linear, curvilinear, or plateau models.

-   The `cldList()` function was revised so that it returns a data frame with the same letter for all groups when there are no significant differences, rather than stop and report "no differences".

-   The `accuracy()` function was revised to include median absolute error.

-   As a side note, all mentions of "sex" as demographic variables in the package were changed to "gender". The same change is being made on <https://rcompanion.org/handbook/> .

## rcompanion v 2.4.6 (2021-11-01)

-   The `efronRSquared()` function was added. It produces Efron's pseudo R-squared from vectors of actual y values, residuals, and predicted values. It can also produce MAE, MAPE, MSE, RMSE, and CV. The `accuracy()` function also produces these statistics directly from a model object.

-   The `countRSquare()` function was added. It produces the count pseudo r-squared for logistic regression or other models with a binary outcome. This statistic is the number of correctly predicted outcomes divided by the total number of outcomes.

-   The `oneSampleDominance()` function was added. It produces a dominance statistic for the one sample case. It can be used as an effect size statistic for the one-sample sign test, and can be used with truly ordinal data.

-   The `pairedSampleDominance()` function was added. This is similar to the `oneSampleDominance()` statistic, but used in the case of two-sample paired data.

-   No longer suggests FSA package.

## rcompanion v 2.4.1 (2021-05-01)

-   The `groupwiseMean()` function was updated to allow for trimmed means, including determining confidence intervals by bootstrap, and to more formally handle NA values.

-   The `nominalSymmetryTest()` function was improved.

-   Documentation was updated for the `groupwiseHuber()` function.

## rcompanion v 2.4.0 (2021-03-23)

-   Two new functions were added that calculate effect sizes for the two-sample case, such as might be used with a Wilcoxon-Mann-Whitney test. These are `wilcoxonPS()`, that calculates Grissom and Kim's probability of superiority, and `wilcoxonOR()`, that calculates Agresti's generalized odds ratio.

-   These functions complement others that could be used as effect size statistics in the two-sample case. `vda()` for Vargha and Delaney's A, `cliffDelta()` for Cliff's delta, `wilcoxonRG()` for Glass rank biserial correlation coefficient, and `wilcoxonR()` for the r statistic Z/sqrt(N).

## rcompanion v 2.3.27 (2021-01-31)

-   Minor bug fix in `freemanTheta()` function.

-   Minor bug fixes, particularly with the `groupwiseCMH()` function to prevent an error in some cases.

## rcompanion v 2.3.25 (2020-02-10)

-   A function was added, `quantileCI()`, that calculates the confidence interval for the median of a numeric or ordered factor variable, based on a method included in Conover, Practical Nonparametric Statistics, 3rd.

## rcompanion v 2.3.21 (2020-01-09)

-   A function was added, `blom()`, which calculates normal scores transformation by Blom, van der Waerden, Tukey, and rankit methods, as well as z-score transformation, and rescaling to a specified range.

-   A function was added, `ordinalEtaSquared()`, that calculates eta-squared as an effect size statistic, following a Kruskal-Wallis test, or for a table with one ordinal variable and one nominal variable; confidence intervals by bootstrap.

-   A function was added, `wilcoxonRG(`), that calculates the Glass rank biserial correlation coefficient effect size for Mann-Whitney two-sample rank-sum test; confidence intervals by bootstrap.

-   A function was added, `wilcoxonPairedRC()`, that calculates the paired-samples rank biserial correlation coefficient effect size for paired Wilcoxon signed-rank test; confidence intervals by bootstrap.

-   A function was added, `wilcoxonOneSampleRC()`, that calculates the paired-samples rank biserial correlation coefficient effect size for one sample Wilcoxon signed-rank test; confidence intervals by bootstrap

-   A modification was made to some effect size functions so that they will not fail when there are instances of the calculation of the statistic failing during the bootstrap procedure. In these cases, the confidence interval is reported as NA (the default) or as the confidence interval based on the successful calculations. (`epsilonSquared()`, `vda()`, `cliffDelta()`, `wilcoxonR()`, `spearmanRho()`, `phi()`, `CramerV()`, `CramerVFit()`, `cohenW()`, `cohenG()`, `freemanTheta()`.

-   A small modification was made to `compareGLM()` and `compareLM()` so that models with a long formula won't cause an error and the function to fail.

-   The function `multiVDA()` was modified to include the Glass rank biserial correlation coefficient.

## rcompanion v 2.3.0 (2019-08-27)

-   The function `phi()` was added to calculate the effect size for a 2 x 2 table of counts, and to optionally produce a confidence interval.

-   The function `cohenW()` was updated to optionally produce a confidence interval.

-   The function `wilcoxonZ()` was added. It extracts the Z statistic from a two-sample, paired, and one-sample Wilcoxon tests.

-   The `wilcoxonR()`, `wilcoxonPairedR()`, and `wilcoxonOneSampleR()` functions were updated to produce confidence intervals more quickly.

## rcompanion v 2.2.2 (2019-07-27)

-   Functions to calculate effect size statistics and confidence intervals were added, specifically for Kendall's w, Spearman's rho, Pearson's r, and Kendall's tau. The relevant functions are: `kendallW()` and `spearmanRho()`.

## rcompanion v 2.2.1 (2019-05-26)

-   Functions to calculate some effect size statistics were updated.

## rcompanion v 2.1.7 (2019-04-09)

-   Confidence intervals by bootstrap for some effect size statistics were added, specifically for Cohen's g, odds ratio, Cohen's P. The relevant function is: `cohenG()`.

-   A function was added to calculate the maximum Vargha and Delaney's A or Cliff's delta as an effect size for Kruskal-Wallis test. `multiVDA()`.

## rcompanion v 2.1.7 (2019-03-02)

-   Confidence intervals by bootstrap for some effect size statistics were added, specifically for Cliff's delta and Vargha and Delaney's A. The relevant functions are: `vda()` and `cliffDelta()`.

## rcompanion v 2.0.10 (2019-01-02)

-   Confidence intervals by bootstrap for some effect size statistics were added, specifically for r for Wilcoxon tests, Cramer's V, epsilon-squared, and Freeman's theta. The relevant functions are: `wilcoxonR()`, `wilcoxonPairedR()`, `wilcoxonOneSampleR()`, `epsilonSquared()`, `freemanTheta()`, `cramerV()`, `cramerVFit()`.

-   The accuracy function was updated to include output for the coefficient of variation for models.

## rcompanion v 2.0.3 (2018-11-01)

-   Some functions were removed that conducted pairwise post-hoc tests where better methods exist: `pairwiseSignTest()`, `pairwiseSignMatrix()`.

-   The `pairwiseMcnemar()` function was updated to include statistics (chi-square, Z) in output.

-   The `transformTukey()` function was modified to allow the user to suppress screen output, and an option was added to return the value of lambda rather than the transformed data.

-   The `nagelkerke()` function was updated to allow the new `lmerModLmerTest` objects

## rcompanion v 2.0.0 (2018-08-13)

-   Some functions were removed that conducted pairwise post-hoc tests where better methods exist. Other functions were modified to decrease the number of imported packages. Removed functions are: `nagelkerkeHermite()`, `pairwiseDifferences()`, `pairwiseOrdinalMatrix()`, `pairwiseOrdinalPairedMatrix()`, `pairwiseOrdinalPairedTest()`, `pairwiseOrdinalTest()`, `pairwiseRobustMatrix()`, `pairwiseRobustTest()`.

## rcompanion v 1.14.0 (2018-05-29)

-   Functions have been amended or added: `multiVDA()`, Calculates Vargha and Delaney's A and Cliff's delta pairwise across groups.

## rcompanion v 1.13.0 (2018-03-25)

Functions have been amended or added:

-   `cohenG()` - Calculates Cohen's g and odds ratio for paired contingency tables, such as those that might be analyzed with McNemar or McNemar-Bowker tests.
-   `cohenH()` - Cohen's h is an effect size to compare two proportions.
-   `wilcoxonR()` - Calculates r effect size for Mann-Whitney, two-sample rank-sum test, or a table with an ordinal variable and a nominal variable with two levels.
-   `wilcoxonOneSampleR()` - Calculates r effect size for a Wilcoxon one-sample signed-rank test.
-   `wilcoxonPairedR()` - Calculates r effect size for a Wilcoxon two-sample paired signed-rank test.
-   `cohenH()` - function to calculate Cohen's h effect size for tables.

## rcompanion v 1.12.3 (2018-02-25)

Functions added:

-   `wilcoxonR()` - function to calculate r effect size for Wilcoxon Mann-Whitney rank sum test

-   `wilcoxonPairedR()` - function to calculate r effect size for Wilcoxon signed rank test

-   `wilcoxonOneSampleR()` - to calculate r effect size for one-sample Wilcoxon signed rank test

Minor improvements to: `freemanTheta()`, `epsilsonSquared()`

## rcompanion v 1.11.3 (2018-02-19)

Improved `cateNelson()` and `cramerV()` functions.

## rcompanion v 1.11.0 (2017-11-12)

-   Added `cramerVFit()` for goodness of fit tests

## rcompanion v 1.10.0 (2017-08-23)

-   Added Efron's pseudo r-square to `accuracy()` function
-   Added `freemanTheta()` function
-   Added `epsilonSquared()` function

## rcompanion v 1.9.1 (2017-07-24)

-   Functions added to test percentiles. `percentileTest()` and `pairwisePercentileTest()` functions.

## rcompanion v 1.8.0 (2017)

-   `scheirerRayHare()` function added

## rcompanion v 1.7.0 (2017)

-   `accuracy()` function added for minimum maximum accuracy, MAPE, RMSE, and NRMSE

-   `cramerV()` function added to calculate Cramer's V for correlation in tables of nominal variables (contingency tables)

## rcompanion v 1.5.1 (2017-04-01)

-   Fixed `nobs` problem in `nagelkerke()` and `nakgelkeHermite()` functions.

## rcompanion v 1.5.0 (2017-01-31)

-   Formula notation added for `cldList()`

## rcompanion v 1.0.0 (2016-09-01)

This package provides custom functions for working through examples and analyses from "Summary and Analysis of Extension Education Program Evaluation in R" and "An R Companion for the Handbook of Biological Statistics".

There are several functions which provide summary statistics for grouped data. These function titles tend to start with `groupwise`". They provide means, medians, geometric means, and Huber M-estimators for groups, along with confidence intervals by traditional methods and bootstrap.

Function titles starting with `pairwise`" conduct pairwise tests among groups as a post-hoc analysis for omnibus tests. At the time of writing, these tests are Mood's median test, sign test (for Friedman test), permutation test, robust anova, and ordinal regression. The output is a table of comparisons and *p*-values, or a matrix of *p*-values that can be parsed into a compact letter display.

There are also functions that are useful for comparing models. `accuracy()`, `compareLM()`, `compareGLM()`, and `pairwiseModelAnova()`. These use goodness-of-fit measures like AIC, BIC, and BICc, or likelihood ratio tests, or accuracy measures like RMSE and Efron's pseudo r-square.

Functions for nominal data include post-hoc tests for Cochran-Mantel-Haenszel test (`groupwiseCMH()`), for McNemar--Bowker (`pairwiseMcnemR()`), and for tests of association like Chi-square, Fisher exact, and G-test (`pairwiseNominalIndependence()`).

A function close to my heart is `cateNelson()`, which performs Cate--Nelson analysis for bivariate data.

Examples and vignettes can be found at <https://rcompanion.org/handbook/> and <https://rcompanion.org/rcompanion/>.
