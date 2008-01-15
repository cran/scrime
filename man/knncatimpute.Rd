\name{knncatimpute}
\alias{knncatimpute}

\title{Missing Value Imputation with kNN}
\description{
  Imputes missing values in a matrix composed of categorical variables
  using \eqn{k} Nearest Neighbors.
}
\usage{
knncatimpute(x, dist = NULL, nn = 3, weights = TRUE)
}

\arguments{
  \item{x}{a numeric matrix containing missing values. All non-missing values
    must be integers between 1 and \eqn{n_{cat}}{n.cat}, where \eqn{n_{cat}}{n.cat}
    is the maximum number of levels the categorical variables in \code{x} can take.
    If the \eqn{k} nearest observations should be used to replace the missing values
    of an observation, then each row must represent one of the observations and each
    column one of the variables. If the \eqn{k} nearest variables should be used
    to impute the missing values of a variable, then each row must correspond to a variable
    and each column to an observation.}
  \item{dist}{either a character string naming the distance measure or a distance matrix.
    If the former, \code{dist} must be either \code{"smc"}, \code{"cohen"}, or \code{"pcc"}.
    If the latter, \code{dist} must be a symmetric matrix having the same number of rows
    as \code{x}. In this case, both the upper and the lower triangle of \code{dist} must
    contain the distances, and the row and column names of \code{dist} must be equal to
    the row names of \code{x}. If \code{NULL}, \code{dist = "smc"} is used.}
  \item{nn}{an integer specifying \eqn{k}, i.e.\ the number of nearest neighbors, used in
    the imputation of the missing values.}
  \item{weights}{should weighted \eqn{k}NN be used to impute the missing values? If \code{TRUE},
    the vote of each nearest neighbor is weighted by the reciprocal of its distance to the observation or variable
    when the missing values of this observation or variable, respectively, are replaced.}
}

\value{
  A matrix of the same size as \code{x} in which all the missing values have been imputed.
}

\references{
Schwender, H.\ (2007). Statistical Analysis of Genotype and Gene Expression Data.
\emph{Dissertation}, Department of Statistics, University of Dortmund.}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{knncatimputeLarge}}, \code{\link{gknn}}, \code{\link{smc}}, \code{\link{pcc}}}

\examples{\dontrun{
# Generate a data set consisting of 200 rows and 50 columns
# in which the values are integers between 1 and 3.
# Afterwards, remove 20 of the values randomly.

mat <- matrix(sample(3, 10000, TRUE), 200)
mat[sample(10000, 20)] <- NA

# Replace the missing values.

mat2 <- knncatimpute(mat)

# Replace the missing values using the 5 nearest neighbors
# and Cohen's Kappa.

mat3 <- knncatimpute(mat, nn = 5, dist = "cohen")

}}

\keyword{classif}
\keyword{NA}
