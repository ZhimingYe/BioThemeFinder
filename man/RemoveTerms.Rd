% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utilities.R
\name{RemoveTerms}
\alias{RemoveTerms}
\title{Remove enriched terms based on regular expression or vector}
\usage{
RemoveTerms(x, Item)
}
\arguments{
\item{x}{a BioThemeFinder object}

\item{Item}{When only one parameter is passed in, the result is filtered using the parameter as a regular expression, and the result matching the regular expression is filtered. When a vector is passed in, all the pathway names contained in the vector are erased from the result.}
}
\value{
A BioThemeFinder object
}
\description{
Remove enriched terms based on regular expression or vector
}
\author{
Zhiming Ye.
}
