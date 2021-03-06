\name{rowMsquares}
\alias{rowMsquares}

\title{Rowwise Linear Trend Test Based on Tables}
\description{
  Given a set of matrices, each representing one group of subjects (e.g.,
  cases and controls in a case-control study), that summarize the numbers 
  of subjects showing the different levels of the categorical variables 
  represented by the rows of the matrices, the value of the linear trend statistic
  based on Pearson's correlation coefficient and described on page 87 of Agresti (2002)
  is computed for each variable.
  
  Using this function instead of \code{rowTrendStats} is in particular recommended
  when the total number of observations is very large.
}
\usage{
rowMsquares(..., listTables = NULL, clScores = NULL, levScores = NULL,
   add.pval = TRUE)
}

\arguments{
  \item{\dots}{numeric matrices in each of which each row corresponds to a ordinal variable
     and each column to one of the ordered levels of these variables. Each of these matrices
     represents one of the groups of interest and comprises the numbers of observations showing 
     the respective levels at the different variables. These matrices can, e.g., generated by
     employing \code{\link{rowTables}}. The dimensions of all matrices must be the same, and 
     the rows and columns must represent the same variables and levels, respectively, in the 
     same order in all matrices. The rowwise sums in a matrix are allowed to differ (which might
     happen if some of the observations are missing for some of the variables.)}
  \item{listTables}{instead of inputting the matrices directly,
     a list consisting of these matrices can be generated and then be used in \code{rowMsquares} 
     by specifying \code{listTables}.}
  \item{clScores}{a numeric vector with one entry for each matrix specifying the score that should
     be assigned to the corresponding group. If \code{NULL}, \code{clScores} is set to \code{1:m},
     where \eqn{m} is the number of groups/matrices, such that the first input matrix (or the first
     entry in \code{listTables}) gets a score of 1, the second a score of 2, and so on.}
  \item{levScores}{a numeric vector with one score for each level of the variables.If not specified, 
     i.e.\ \code{NULL}, the column names of the matrices are interpreted as scores.}
  \item{add.pval}{should p-values be added to the output? If \code{FALSE}, only the
     rowwise values of the linear trend test statistic will be returned. If \code{TRUE},
     additionally the (raw) p-values based on an approximation to the ChiSquare-distribution
     with 1 degree of freedom are returned.}
}
\details{
 This is an extension of the Cochran-Armitage trend test from two to several classes. The
 statistic of the Cochran-Armitage trend test can be obtained by multiplying the statistic
 of this general linear trend test with \eqn{n / (n - 1)}, where \eqn{n} is the number
 of observations.
}
\value{
  Either a vector containing the rowwise values of the linear trend test statistic
  (if \code{add.pval = FALSE}), or a list containing these values (\code{stats}),
  and the (raw) p-values (\code{rawp}) not adjusted for multiple comparisons (if \code{add.pval = TRUE}).
}

\references{ 
  Agresti, A.\ (2002). \emph{Categorical Data Analysis}. Wiley, Hoboken, NJ. 2nd Edition.

  Mantel, N.\ (1963). Chi-Square Test with one Degree of Freedom: Extensions of the Mantel-Haenszel
  Procedure. \emph{Journal of the American Statistical Association}, 58, 690-700.
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}
\note{
  The usual contingency table for a variable can be obtained from the matrices
  by forming a variable-specific matrix in which each row consists of the
  row of one of these matrices.
}
\seealso{\code{\link{rowTrendStats}}, \code{\link{rowCATTs}}, \code{\link{rowChisqMultiClass}}}
\examples{\dontrun{
# Generate a matrix containing data for 10 categorical 
# variables with levels 1, 2, 3.

mat <- matrix(sample(3, 500, TRUE), 10)

# Now assume that we consider a case-control study,
# i.e. two groups, and that the first 25 columns 
# of mat correspond to cases and the remaining 25 
# columns to cases. Then a vector containing the 
# class labels is given by

cl <- rep(1:2, e=25)

# and the matrices summarizing the numbers of subjects
# showing the respective levels at the different variables
# are computed by

cases <- rowTables(mat[, cl==1])
controls <- rowTables(mat[,cl==2])

# The values of the rowwise liner trend test are
# computed by

rowMsquares(cases, controls)

# which leads to the same results as

listTab <- list(cases, controls)
rowMsquares(listTables = listTab)

# or as

rowTrendStats(mat, cl, use.n = FALSE)

# or as

out <- rowCATTs(cases, controls)
n <- ncol(mat)
out$stats * (n - 1) / n

 }
}

\keyword{htest}
\keyword{array}
