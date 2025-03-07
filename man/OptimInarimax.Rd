% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{OptimInarimax}
\alias{OptimInarimax}
\title{A improved parameters estimation method for TSA::arimax}
\usage{
OptimInarimax(y, order, xtransf, transfer)
}
\arguments{
\item{y}{a univariate time series}

\item{order}{regular ARIMA order}

\item{xtransf}{intervention variables}

\item{transfer}{a list consisting of the ARMA orders for each transfer (distributed lag) covariate.}
}
\value{
The ArimaOptExtr function returns a fitted model with robust parameter estimates.
}
\description{
A improved parameters estimation method for TSA::arimax
}
