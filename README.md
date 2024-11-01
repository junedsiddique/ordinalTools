
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ordinalTools: Tools for working with output from ordinal regression models

## Description

<strong>ordinalTools</strong> is a package for working with output from
ordinal regression models. As currently implemented, the package uses
output from ordinal logistic regressions fit using the \`polr()’
function from the package MASS to calculate means based on
user-specified covariate constrasts. Standard errors are estimated
either using the delta method or bootstrapping.

## Basic Features

- The package contains a single function named `ordinal_means()` where
  the only required argument is `object`, an object of class polr.

- Methods for standard generics are provided, i.e., `coef()`,
  `summary()`, `confint()`, and `print()`.

- Bootstrap standard errors and confidence intervals can be estimated
  using se.type=“bootstrap” and by specifying the data frame used to fit
  the original polr regression model.

- If two contrasts are specified, the function will also calculate the
  difference in means based on the two contrasts as well as its standard
  error and confidence interval.

## Basic Use

Let `y` denote an ordinal outcome, and `x1` and `x2` covariates. An
ordinal logistic regression model model with `y` as outcome and `x1` and
`x2` as covariates is fitted using the polr package

``` r
library(MASS)
polr.fit <- polr(y ~ x1 + x2, data = DF, Hess=TRUE)
```

Using the polr.fit object, we can calculate the mean of `y` for various
levels of the covariates `x1` and `x2`

``` r
fit <- ordinal_means(polr.fit, contrast1=c(0,0), contrast2=c(1,0))

summary(fit)
confint(fit)
```

To calculate bootstrap standard errors and bootstrap confidence
intervals, one must specify the data frame used in the polr call.

``` r
fit <- ordinal_means(polr.fit, data=DF, contrast1=c(0,0), contrast2=c(1,0),
                     se.type="bootstrap", R=100, conf.level = 0.95)

summary(fit)
confint(fit)
```

## Installation

The development version of the package can be installed from GitHub
using the **devtools** package:

``` r
devtools::install_github("junedsiddique/ordinalTools")
```
