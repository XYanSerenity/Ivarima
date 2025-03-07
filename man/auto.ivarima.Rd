% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{auto.ivarima}
\alias{auto.ivarima}
\title{An automated, rule-based approach for identifying the optimal intervention impact model in interrupted time series,
constrained to cases where the intervention variable is an impulse function, the impact model follows
a first-order intervention transfer function, and the noise component is a non-seasonal ARIMA model}
\usage{
auto.ivarima(y, its_start, LTFmaxk)
}
\arguments{
\item{y}{is a time series.}

\item{its_start}{a number between 1 and length(y)-1 stating the time point of the start of the intervention.}

\item{LTFmaxk}{The LTFmaxk accounts for the longest time-lagged response in Linear Transfer Function  model.}
}
\value{
The function returns the optimal model fitting results found.
}
\description{
An automated, rule-based approach for identifying the optimal intervention impact model in interrupted time series,
constrained to cases where the intervention variable is an impulse function, the impact model follows
a first-order intervention transfer function, and the noise component is a non-seasonal ARIMA model
}
\examples{
w0 <- 0.5
n <- 100
xt6 <- c(rep(0, 56), rep(1, 44))
e <- arima.sim(model = list(order = c(1, 0, 0), ar = 0.4, ma = NULL), n = 100, sd = 0.1)
xxt <- w0 * xt6
mu <- 4 + xxt
y <- mu + e
model <- auto.ivarima(y = y, its_start = 51, LTFmaxk = 10)
summary(model)
}
