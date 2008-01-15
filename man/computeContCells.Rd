\name{computeContCells}
\alias{computeContCells}

\title{Pairwise Contingency Tables}
\description{
  Computes a contingency table for each pair of rows of a matrix, and stores
  all contigency table in a matrix.
}
\usage{
computeContCells(data, computeExp = TRUE, justDiag = FALSE, 
   check = TRUE, n.cat = NULL)
}

\arguments{
  \item{data}{a numeric matrix consisting of integers between 1 and \code{n.cat}. 
    It is assumed that each row of these matrix represents
    a variable. Missing values and different numbers of categories a variable can take are
    allowed.}
  \item{computeExp}{should the numbers of observations expected under the null hypothesis that
    the respective two variables are independent also be computed? Required
    when \code{computeContCells} is used to compute Pearson's \eqn{\chi^2}{ChiSquare}-statistic.}
  \item{justDiag}{should only the diagonal elements of the contingency tables,
    i.e.\ \eqn{n_{ii}}{n.ii}, \eqn{i = 1, \dots}, \code{n.cat}, be computed?}
  \item{check}{should \code{data} be checked more thoroughly? It is highly recommended to use
    \code{check = TRUE}.}
  \item{n.cat}{integer specifying the maximum number of levels a variable can take. If \code{NULL},
    this number will be computed. It is highly recommended not to change the default.}
}

\value{
  A list consisting of two matrices each consisting of \eqn{m * (m - 1) / 2} rows and 
  \code{n.cat}\eqn{^2}{^2} columns, where \eqn{m} is the number of rows of \code{data}. 
  One of these matrices
  called \code{mat.obs} contains in each row the values of the contingency table for
  a particular pair of rows of \code{data}, where the contigency table of the variables represented
  by the \emph{i}th and \emph{j}th row of \code{data} is shown in the
  \eqn{j + m * (i - 1) - i * (i - 1) / 2} row of \code{mat.obs}.
  The other matrix called \code{mat.exp} consists of
  the corresponding numbers of observations expected under the null hypothesis that
  the respective two variables are independent.
}

\references{
Schwender, H.\ (2007). A Note on the Simultaneous Computation of Thousands of 
Pearson's \eqn{\chi^2}{ChiSquare}-Statistics. \emph{Technical Report}, SFB 475,
Deparment of Statistics, University of Dortmund. 
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{computeContClass}}, \code{\link{rowChisqStats}}}

\examples{\dontrun{
# Generate an example data set consisting of 5 rows (variables)
# and 200 columns (observations) by randomly drawing integers 
# between 1 and 3.

mat <- matrix(sample(3, 1000, TRUE), 5)

# Generate the matrix containing the contingency tables for each
# pair of rows of mat.

out <- computeContCells(mat)

# out contains both the observed numbers of observations 
# summarized by contingency tables

out$mat.obs

# and the number of observations expected under the null hypothesis
# of independence.

out$mat.exp

# If, e.g., only the observed number of observations having the same 
# value is of interest, call

computeContCells(mat, computeExp = FALSE, justDiag = TRUE)

}}

\keyword{manip}
