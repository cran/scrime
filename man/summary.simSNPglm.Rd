\name{summary.simSNPglm}
\alias{summary.simSNPglm}

\title{Summarizing a simSNPglm object}
\description{
  Summarizes an object of class \code{simSNPglm}.
}
\usage{
\method{summary}{simSNPglm}(object, digits = 3, ...)
}

\arguments{
  \item{object}{an object of class \code{simSNPglm}, i.e.\ the output of \code{\link{simulateSNPglm}}.}
  \item{digits}{number of digits used in the output.}
  \item{\dots}{Ignored.}
}
\value{
  Shows the model used in \code{\link{simulateSNPglm}} to generate the values of the response.
  If the response is binary, then it additionally shows and returns a contingency table of the
  numbers of cases and controls and the probabilities for being a case.
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{simulateSNPglm}}}

\examples{\dontrun{
# The default simulated data set.

sim1 <- simulateSNPglm()
sim1

# A bit more information: Table of probability of being a case
# vs. number of cases and controls.

summary(sim1)

}}

\keyword{print}
