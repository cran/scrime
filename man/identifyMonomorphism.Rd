\name{smc}
\alias{identifyMonomorphism}

\title{Identification of Constant Variables}
\description{
  Identifies the rows of a matrix that only show one level.
}
\usage{
identifyMonomorphism(x)
}

\arguments{
  \item{x}{a matrix consisting of integers between 1 and \eqn{n_{cat}}{n.cat},
    where \eqn{n_{cat}}{n.cat} is the number of levels the variables in \code{x}
    can take. Each row of \code{x} is assumed to represent one variable. Missing values are allowed.}
}

\value{
  Numeric vector containing the rows of \code{x} showing only one level.
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\examples{\dontrun{
# Generate a data set consisting of 10 rows and 15 columns,
# where the values are randomly drawn from the integers 1, 2, and 3,
# and row 3 and 7 consist only of one level.

mat <- matrix(sample(3, 2000, TRUE), 10)
mat[3, ] <- 1
mat[7, ] <- 2

identifyMonomorphism(mat)

}}

\keyword{manip}
