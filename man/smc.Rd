\name{smc}
\alias{smc}
\alias{cohen}

\title{Simple Matching Coefficient and Cohen's Kappa}
\description{
  Computes the values of (or the distance based on) the simple matching coefficient
  or Cohen's Kappa, respectively, for each pair of rows of a matrix.
}
\usage{
smc(x, dist = FALSE)
cohen(x, dist = FALSE)
}

\arguments{
  \item{x}{a matrix consisting of integers between 1 and \eqn{n_{cat}}{n.cat},
    where \eqn{n_{cat}}{n.cat} is the number of levels the variables in \code{x}
    can take. Missing values are allowed.}
  \item{dist}{should the distance based on the simple matching coefficient or Cohen's Kappa, respectively,
    be computed? Note that, e.g., \code{smc(x, dist = TRUE)} is equal to \code{1 - smc(x, dist = FALSE)}.}
}

\value{
  A matrix with \code{nrow(x)} columns and rows containing the distances or similarities.
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{pcc}}}

\examples{\dontrun{
# Generate a data set consisting of 10 rows and 200 columns,
# where the values are randomly drawn from the integers 1, 2, and 3.

mat <- matrix(sample(3, 2000, TRUE), 10)

# For each pair of row, the value of the simple matching coefficient
# can be obtained by

smc(mat)

# and the distance based on the SMC by

smc(mat, dist = TRUE)

}}

\keyword{manip}
