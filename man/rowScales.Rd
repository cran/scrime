\name{rowScales}
\alias{rowScales}

\title{Rowwise Scaling}
\description{
  Scales each row of a matrix such that the values in this row have
  zero mean and a standard deviation of 1.
}
\usage{
rowScales(X, add.stats = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric matrix whose rows should be scaled. Missing values are allowed.}
  \item{add.stats}{should the rowwise means and standard deviations of \code{X}
    be returned?}
}

\value{
  If \code{add.stats = FALSE}, a matrix of the same dimensions as \code{X}
  with a rowwise mean of zero and a rowwise standard deviation of 1.
  If \code{add.stats = TRUE}, a list containing this matrix and the
  rowwise means and standard deviations of the input matrix.
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{rowCors}}}
\examples{\dontrun{
# Generate a matrix containing data for 10 categorical 
# variables with levels 1, 2, 3.

mat <- matrix(sample(3, 500, TRUE), 10)

rowScales(mat)

}}
\keyword{array}
\keyword{manip}
