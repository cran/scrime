\name{rowCors}
\alias{rowCors}

\title{Rowwise Correlation with a Vector}
\description{
  Computes Pearson's correlation coefficient of a vector with each row of a matrix.
}
\usage{
rowCors(X, y, trendStat = FALSE, use.n = NULL)
}

\arguments{
  \item{X}{a numeric matrix in which each row represents a variable and each
     column an observation.}
  \item{y}{a numeric vector of length \code{ncol(X)}.}
  \item{trendStat}{instead of the correlation coefficients should the values of
    the statistic for a test of linear trend based on this coefficient be returned?
    If \code{TRUE}, then it is assumed that all variables in \code{X} and the 
    variable represented by \code{y} are ordinal, and the values in \code{X} and \code{y}
    represent scores for the different levels.}
  \item{use.n}{should the squared values of the correlation coefficient be multiplied
    by \code{ncol(X)}? Ignored if \code{trendStat = FALSE}. If \code{FALSE}, the
    squared values are multiplied by \code{ncol(X) - 1}. By default, the squared values
    are multiplied by \code{ncol(X)} if \code{y} shows two levels, leading to
    the Cochran-Armitage test of trend. Otherwise, they are multiplied by \code{ncol(X) - 1}.}
}

\value{
  A vector containing the rowwise values of Pearson's correlation coefficient (if
  \code{trendStat = FALSE} or the rowwise values of the trend statistics (if
  \code{trendStat = TRUE}.
}
\references{ 
  Agresti, A.\ (2002). \emph{Categorical Data Analysis}. Wiley, Hoboken, NJ. 2nd Edition.
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{rowTrendStats}}, \code{\link{rowCATTs}}, \code{\link{rowMsquares}}}

\examples{\dontrun{
# Generate a random matrix containing 10 continuous variables
# and a vector representing a continuous variable.

mat <- matrix(runif(200, 0, 20), 10)
y <- sample(runif(20, 0, 20))

# The correlations between y and each of row of mat are
# computed by

rowCors(mat, y)

# Generate a random binary vector and a matrix consisting
# of 10 ordinal variables with levels 0, 1, 2, where these 
# values can be interpreted as scores for the differ 
# categories.

mat <- matrix(sample(0:2, 500, TRUE), 10)
y <- sample(0:1, 50, TRUE)

# The values of the Cochran-Armitage trend statistic are
# computed by

rowCors(mat, y, trendStat = TRUE)

# If the values of the general test of linear trend described
# on page 87 of Agresti (2002) should be computed, then call

rowCors(mat, y, trendStat = TRUE, use.n = FALSE) 

}}

\keyword{manip}
\keyword{array}
