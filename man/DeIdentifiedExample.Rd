\name{DeIdentifiedExample}
\alias{DeIdentifiedExample}
\docType{data}
\title{
A simple example for testing pairwise comparison methods
}
\description{
A simple example of data for testing pairwise comparison methods composed of a normally distributed response and a factor with 7 levels.
}
\usage{data("DeIdentifiedExample")}
\format{
  A data frame with 70 observations on the following 2 variables.
  \describe{
    \item{\code{Y}}{a numeric vector containing a response variable}
    \item{\code{Treatment}}{a factor with 7 levels \code{A} \code{B} \code{C} \code{D} \code{E} \code{F} \code{G}}
  }
}
\details{
 De identified real dataset for testing
}
\source{
A real data set with altered names, so that the source (irrelevant in this context) cannot be identified.
}
\examples{
data(DeIdentifiedExample)
str(DeIdentifiedExample)
}
\keyword{datasets}
