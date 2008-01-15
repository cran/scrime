\name{computeContClass}
\alias{computeContClass}

\title{Rowwise Contigency Tables}
\description{
  Generates a matrix containing a contingency table for each row of a matrix
  and a vector of class labels.
}
\usage{
computeContClass(data, cl, n.cat)
}

\arguments{
  \item{data}{a numeric matrix consisting of integers between 1 and \code{n.cat}. Each
    row of data is assumed to represent a variable, and each column to represent an
    observation. Missing values are not allowed. All variables must
    comprise the same number of levels.}
  \item{cl}{a numeric vector of length \code{ncol(data)} specifying the
    class labels of the observations represented by the columns of data. \code{cl} must
    consist of integers between 1 and \eqn{n_{cl}}{n.cl}, where \eqn{n_{cl}}{n.cl} is
    the number of classes.}
  \item{n.cat}{an integer giving the number of levels the variables can take. If not
    specified, \code{n.cat} will be determined automatically. It is highly
    recommended not to specify \code{n.cat}.}
}

\value{
  A list composed of the following two matrices:
  \item{mat.obs}{a matrix consisting of \eqn{m} rows and \code{n.cat} \eqn{* n_{cl}}{* n.cl} columns,
    where \eqn{m} is the number of variables, i.e.\ the number of rows of \code{data}. 
    Each row of data shows the contingency table of \code{cl} and the corresponding row of \code{data}.}
  \item{mat.exp}{a matrix of the same size as \code{mat.obs} containing the numbers of observations
    expected under the null hypothesis of equal distribution in all classes that correspond to the
    respective entries in \code{mat.obs}.}
}

\references{
Schwender, H.\ (2007). A Note on the Simultaneous Computation of Thousands of 
Pearson's \eqn{\chi^2}{ChiSquare}-Statistics. \emph{Technical Report}, SFB 475,
Deparment of Statistics, University of Dortmund. 
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{computeContCells}}, \code{\link{rowChisqStats}}}

\examples{\dontrun{
# Generate an example data set consisting of 10 rows (variables)
# and 200 columns (observations) by randomly drawing integers 
# between 1 and 3, and a vector of class labels of length 200
# indicating that the first 100 observation belong to class 1
# and the other 100 to class 2. 

mat <- matrix(sample(3, 2000, TRUE), 10)
cl <- rep(1:2, e = 100)

# Applying computeContClass to this data set

out <- computeContClass(mat, cl)

# generates the observed numbers of observations

out$mat.obs

# and the corresponding expected numbers of observations.

out$mat.exp

}}


\keyword{manip}
