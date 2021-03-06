% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GroupClustering.R
\name{MedianBootCI}
\alias{MedianBootCI}
\title{Confidence intervals for medians}
\usage{
MedianBootCI(y, g, Nboots = 9999, Coverage = 0.95, digits = 4)
}
\arguments{
\item{y}{is a vector with the response.}

\item{g}{is a single classification factor.}

\item{Nboots}{number of bootstrap samples (default = 9999).}

\item{Coverage}{the coverage probability of the confidence interval.}

\item{digits}{the number of digits (devfault = 4)
(default = 0.95)
"none" does not adjust}
}
\value{
a matrix with three columns, the median, the lower bound and the
upper bound of the confidence interval and one row for each level of g.
}
\description{
Confidence intervals for medians
}
\examples{
MedianBootCI (y = DeIdentifiedExample$Y,
                                  g = DeIdentifiedExample$Treatment)
}
\author{
Rodrigo Labouriau
}
