% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MGM.r
\name{MGM}
\alias{MGM}
\title{Defining S3 object "MGM"}
\usage{
MGM(X, Y, g)
}
\arguments{
\item{X}{data frame or matrix of continuous variables (row: observation, column: variable)}

\item{Y}{data frame or matrix of discrete variables (row: observation, column: variable)}

\item{g}{group index, needed for temporary files}
}
\value{
An S3 `MGM` object, containing data, network parameters, and the 1st derivatives
}
\description{
Defining S3 object "MGM"
}
