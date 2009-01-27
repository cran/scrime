\name{rowTrendStats}
\alias{rowTrendStats}

\title{Rowwise Linear Trend Tests}
\description{
  Computes for each row of a matrix the value of the statistic of a linear trend test
  for testing whether the ordinal variable corresponding to the row of the matrix
  is associated with an ordinal response. 
  
  In the two-class case, the statistic of the Cochran-Armitage trend test is computed
  by default.
}
\usage{
rowTrendStats(X, y, use.n = NULL, add.pval = TRUE)
}

\arguments{
  \item{X}{a numeric matrix in which each row represents an ordinal variable and
    each column corresponds to an observation. The entries of this matrix are
    interpreted as scores for the different (ordered) levels of the respective variables.}
  \item{y}{a numeric vector of length \code{ncol(X)} containing the class labels
    of the observations represented by the columns of \code{X}, where these labels
    are interpreted as scores for the different classes.}
  \item{use.n}{should the squared values of Pearson's correlation coefficient be multiplied
    by \code{ncol(X)} to generate the values of the test statistic? If \code{FALSE}, the
    squared values are multiplied by \code{ncol(X) - 1}. By default, the squared values
    are multiplied by \code{ncol(X)} if \code{y} shows two levels, leading to
    the Cochran-Armitage test of trend. Otherwise, they are multiplied by \code{ncol(X) - 1}
    leading to the linear trend test statistic of Mantel (1963) described, e.g., on page
    87 of Agresti (2002).}
  \item{add.pval}{should p-values be added to the output? If \code{FALSE}, only the
     rowwise values of the linear trend test statistic will be returned. If \code{TRUE},
     additionally the (raw) p-values based on an approximation to the ChiSquare-distribution
     with 1 degree of freedom are returned.}
}

\references{ 
  Agresti, A.\ (2002). \emph{Categorical Data Analysis}. Wiley, Hoboken, NJ. 2nd Edition.

  Armitage, P.\ (1955). Tests for Linear Trends in Proportions and Frequencies. 
  \emph{Biometrics}, 11, 375-386.
  
  Cochran, W.~G.\ (1954). Some Methods for Strengthening the Common ChiSquare Tests.
  \emph{Biometrics}, 10, 417-451.
  
  Mantel, N.\ (1963). Chi-Square Test with one Degree of Freedom: Extensions of the Mantel-Haenszel
  Procedure. \emph{Journal of the American Statistical Association}, 58, 690-700.
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}
\seealso{\code{\link{rowMsquares}}, \code{\link{rowCATTs}}, \code{\link{rowChisqMultiClass}}}
\examples{\dontrun{
# Generate a matrix containing data for 10 categorical 
# variables with levels 1, 2, 3.

mat <- matrix(sample(3, 500, TRUE), 10)

# Now assume that the first 25 columns correspond to
# cases and the remaining 25 columns to cases. Then
# a vector containing the class labels is given by

cl <- rep(0:1, e=25)

# The values of the Cochran-Armitage trend test can
# then be computed by

rowTrendStats(mat, cl)

# This leads to the same results as

cases <- rowTables(mat[, cl==1])
controls <- rowTables(mat[,cl==0])

rowCATTs(cases, controls)

# or as

out <- rowMsquares(cases, controls)
n <- ncol(mat)
out$stats * n / (n - 1)

  }
}

\keyword{htest}
\keyword{array}
