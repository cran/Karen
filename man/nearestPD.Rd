% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{nearestPD}
\alias{nearestPD}
\title{Nearest Positive Definite Matrix}
\usage{
nearestPD(A, ...)
}
\arguments{
\item{A}{numeric \eqn{n \times n}{n by n} approximately positive definite matrix, typically an approximation to a correlation or covariance matrix.
If A is not symmetric (and ensureSymmetry is not false), symmpart(A) is used.}

\item{...}{Further arguments to be passed to nearPD (see package Matrix for details).}
}
\value{
The nearest positive definite matrix of A.
}
\description{
This function first check if a matrix A is positive definite, typically a correlation or variance-covariance matrix.
If A is not positive definite, this function computes the nearest positive definite matrix of A using the function nearPD from package Matrix.
}
\examples{
nearestPD(diag(c(1,0,1)))
}
