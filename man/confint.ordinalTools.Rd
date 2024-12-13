% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{confint.ordinalTools}
\alias{confint.ordinalTools}
\title{Confidence Intervals for Means}
\usage{
\method{confint}{ordinalTools}(object, ...)
}
\arguments{
\item{object}{a object of class "ordinalTools"}

\item{...}{other arguments}
}
\value{
A matrix (or vector) with columns giving
lower and upper confidence limits for each parameter.
These will be labelled as (1-level)/2 and 1 - (1-level)/2 in \%
}
\description{
Confidence Intervals for Means
}
\examples{
\dontrun{
fit <- ordinal_means(polr.fit)
confint(fit)
}
}
