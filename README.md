
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Ivarima

<!-- badges: start -->
<!-- badges: end -->

This package implements a rule-based automated method for identifying
the optimal intervention impact model for intervention ARIMA model.This
method utilizes the impulse response function of the LTF model to infer
the possible delay time and adopts the sequential Koyck model
identification rules to guide the identification of decay pattern (h)
and unpatterned spikes (r). An example describing how to install and use
this package is described below. A more detailed tutorial, including the
data analysis described in the paper, is also available with this
package (pdf file)

## Installation

You can install the development version of Ivarima from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("XYanSerenity/Ivarima")
```

  The ‘Ivarima’ R package includes four main functions:
  1. IdenPatternb() – Identify the candidate b values. This function determines the b values that best align with the theoretical b-pattern.
  2. OptimInarimax() – An improved parameter estimation method for TSA::arimax. This function focuses on obtaining suitable initial parameter estimates for the intervention component and intercept term, which serve as starting values for the BFGS algorithm search, reducing the risk of suboptimal results.
  3. FindMinAic() – Selects the optimal model from the evaluated interventional ARIMA models using the evaluation criteria we developed.
  4. auto.ivarima() – Automates the identification process of the impact pattern.


## Example

This is a basic example which shows you how to solve a common problem:

``` r
#library(Ivarima)

library(Ivarima)
w0 <- 0.5
n <- 100
xt6 <- c(rep(0, 56), rep(1, 44))
e <- arima.sim(model = list(order = c(1, 0, 0), ar = 0.4, ma = NULL), n = 100, sd = 0.1)
xxt <- w0 * xt6
mu <- 4 + xxt
y <- mu + e
model <- auto.ivarima(y = y, its_start = 51, LTFmaxk = 10)
summary(model)
```

The auto.ivarima function returns a list with two elements: (i) identified impact pattern, a character string where the first three digits represent the orders of h, b, and r, respectively, the remaining characters “step” represent a step intervention, “pulse” represent a pulse intervention, “pulseD” represent a pulse with a decay intervention;  (ii) the fitted optimal interventional ARIMA model.
