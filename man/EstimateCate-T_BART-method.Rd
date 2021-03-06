% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TBART.R
\name{EstimateCate,T_BART-method}
\alias{EstimateCate,T_BART-method}
\title{EstimateCate-T_BART-method}
\usage{
\S4method{EstimateCate}{T_BART}(theObject, feature_new, verbose = FALSE, return_CI = FALSE)
}
\arguments{
\item{theObject}{A `MetaLearner` object.}

\item{feature_new}{A feature data frame.}

\item{verbose}{TRUE for detailed output FALSE for no output}

\item{return_CI}{If TRUE, return predictions and their confidence intervals;}
}
\description{
Returns the estimated CATE.
}
