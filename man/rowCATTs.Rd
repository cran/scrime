\name{rowCATTs}
\alias{rowCATTs}

\title{Rowwise Cochran-Armitage Trend Test Based on Tables}
\description{
  Given two matrices, each representing one group of subjects (e.g.,
  cases and controls in a case-control study), that summarize the numbers 
  of subjects showing the different (ordered) levels of the ordinal variables 
  represented by the rows of the matrices, the value of the Cochran-Armitage
  Trend Test statistic is computed for each variable.
  
  Using this function instead of \code{rowTrendStats} is in particular recommended
  when the total number of observations is very large.
}
\usage{
rowCATTs(cases, controls, scores = NULL, add.pval = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cases}{a numeric matrix in which each row represents one ordinal
     variable and each column one of the ordered levels that the variables exhibit. The
     entries of this matrix are the numbers of observations from one group (e.g.,
     the cases in a case-control study) showing a particular
     level at the different variables. Such a matrix can, e.g., be generated
     by \code{\link{rowTables}}. The rowwise sums of \code{cases} are allowed to
     differ between variables (which might happen when some of the observations
     are missing for some of the variables).}
  \item{controls}{a numeric matrix of the same dimensions as \code{cases} comprising
     the numbers of observations from the second group (e.g., the controls in a
     case-control study) that show the respective level at the different ordinal variables.
     The rows of \code{controls} must represent the same variables in the same
     order as \code{cases}, and the columns must represent the same levels in the
     same order. This matrix can also be generated by employing \code{\link{rowTables}}.
     The rowwise sums of \code{controls} are allowed to
     differ between variables (which might happen when some of the observations
     are missing for some of the variables).}
  \item{scores}{a numeric vector of length \code{ncol(cases)} containing the scores
     for the different levels. If not specified, i.e.\ \code{NULL}, the column names
     of \code{cases} are interpreted as scores.}
  \item{add.pval}{should p-values be added to the output? If \code{FALSE}, only the
     rowwise values of the Cochran-Armitage trend test statistic will be returned. If \code{TRUE},
     additionally the (raw) p-values based on an approximation to the ChiSquare-distribution
     with 1 degree of freedom are returned.}
}

\value{
  Either a vector containing the rowwise values of the Cochran-Armitage trend test statistic
  (if \code{add.pval = FALSE}), or a list containing these values (\code{stats}),
  and the (raw) p-values (\code{rawp}) not adjusted for multiple comparisons (if \code{add.pval = TRUE}).
}
\references{ 
  Agresti, A.\ (2002). \emph{Categorical Data Analysis}. Wiley, Hoboken, NJ. 2nd Edition.

  Armitage, P.\ (1955). Tests for Linear Trends in Proportions and Frequencies. 
  \emph{Biometrics}, 11, 375-386.
  
  Cochran, W.~G.\ (1954). Some Methods for Strengthening the Common ChiSquare Tests.
  \emph{Biometrics}, 10, 417-451.
}
\author{Holger Schwender, \email{holger.schwender@udo.edu}}
\note{
  The usual contingency table for a variable can be obtained from the matrices
  by forming a variable-specific matrix in which each row consists of the
  row of one of these matrices.
}
\seealso{\code{\link{rowTrendStats}}, \code{\link{rowMsquares}}, \code{\link{rowChisq2Class}}}
\examples{\dontrun{
# Generate a matrix containing data for 10 categorical 
# variables with levels 1, 2, 3.

mat <- matrix(sample(3, 500, TRUE), 10)

# Now assume that the first 25 columns correspond to
# cases and the remaining 25 columns to cases. Then
# a vector containing the class labels is given by

cl <- rep(1:2, e=25)

# and the matrices summarizing the numbers of subjects
# showing the respective levels at the different variables
# are computed by

cases <- rowTables(mat[, cl==1])
controls <- rowTables(mat[,cl==2])

# The values of the rowwise Cochran-Armitage trend test 
# are computed by

rowCATTs(cases, controls)

# which leads to the same results as

rowTrendStats(mat, cl)

# or as 

out <- rowMsquares(cases, controls)
n <- ncol(mat)
out$stats * n / (n-1)

 }
}

\keyword{htest}
\keyword{array}
