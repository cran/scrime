\name{pcc}
\alias{pcc}

\title{Pearson's Contingency Coefficient}
\description{
  Computes the values of (the corrected) Pearson's contingency coefficient for
  all pairs of rows of a matrix. 
}
\usage{
pcc(x, dist = FALSE, corrected = TRUE, version = 1)
}

\arguments{
  \item{x}{a numeric matrix consisting of integers between 1 and \eqn{n_{cat}}{n.cat},
    where \eqn{n_{cat}}{n.cat} is the maximum number of levels a variable in \code{x} can take.}
  \item{dist}{should the distance based on Pearson's contingency coefficient be computed?
    For how this distance is computed, see \code{version}.}
  \item{corrected}{should Pearson's contingency coefficient be corrected such that it can
    take values between 0 and 1? If not corrected, it takes values between and 0
    and \eqn{\sqrt{(a - 1) / a}}{sqrt((a - 1) / a)}, 
    where \eqn{a} is the minimum of the numbers of levels that the respective
    two variables can take. Must be set to \code{TRUE}, if \code{dist = TRUE}.}
  \item{version}{a numeric value -- either 1, 2, or 3 -- specifying how the distance is computed.
    Ignored if \code{dist = FALSE}. If \code{1}, \eqn{\sqrt{1 - Cont^2}}{sqrt(1 - Cont^2)} is computed,
    where \eqn{Cont} denotes Pearson's contigency coefficient. If 2, \eqn{1 - Cont} is determined, and if 3,
    \eqn{1 - Cont^2}{1 - Cont^2} is returned.}
}

\value{
  A matrix with \code{nrow(x)} columns and rows containing the values of (or distances based on)
  the (corrected) Pearson's contigency coefficient for all pairs of rows of \code{x}.
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{smc}}}

\examples{\dontrun{
# Generate a data set consisting of 10 rows and 200 columns,
# where the values are randomly drawn from the integers 1, 2, and 3.

mat <- matrix(sample(3, 2000, TRUE), 10)

# For each pair of rows of mat, the value of the corrected Pearson's 
# contingency coefficient is then obtained by

out1 <- pcc(mat)
out1

# and the distances based on this coefficient by

out2 <- pcc(mat, dist = TRUE)
out2

# Note that if version is set to 1 (default) in pcc, then

all.equal(sqrt(1 - out1^2), out2)

}}



\keyword{manip}
