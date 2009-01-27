\name{rowChisqStats}
\alias{rowChisqStats}

\title{Rowwise Pearson's ChiSquare Statistic}
\description{
  Computes for each row of a matrix the value of Pearson's ChiSquare statistic
  for testing if the corresponding categorical variable is associated with a (categorical)
  response, or determines for each pair of rows of a matrix the value of
  Pearson's ChiSquare statistic for testing if the two corresponding variables
  are independent.
}
\usage{
rowChisqStats(data, cl, compPval = TRUE, asMatrix = TRUE)
}

\arguments{
  \item{data}{a numeric matrix consisting of the integers between 1 and \eqn{n_{cat}}{n.cat},
    where \eqn{n_{cat}}{n.cat} is the maximum number of levels the categorical variables can
    take. Each row of \code{data} must correspond to a variable, each row to an observation.
    Missing values and different numbers of levels a variable might take are allowed.}
  \item{cl}{a numeric vector of length \code{ncol(data)} containing the class labels for
   the observations represented by the columns of \code{data}. The class labels must be
   coded by the integers between 1 and \eqn{n_{cl}}{n.cl}, where \eqn{n_{cl}}{n.cl}
   is the number of classes. If missing, the value of the statistic for Pearson's
   \eqn{\chi^2}{ChiSquare}-test of independence will be computed for each pair of rows
   of \code{data}. Otherwise, the value of Pearson's \eqn{\chi^2}{ChiSquare}-statistic
   for testing if the distribution of the variable differs between the groups specified
   by \code{cl} will be determined for each row of \code{data}.}
  \item{compPval}{should also the p-value (based on the approximation to a 
    \eqn{\chi^2}{ChiSquare}-distribution) be computed?}
  \item{asMatrix}{should the pairwise test scores be returned as matrix? Ignored
    if \code{cl} is specified. If \code{TRUE},
    a matrix with \eqn{m} rows and columns is returned that contains the values
    of Pearson's \eqn{\chi^2}{ChiSquare}-statistic in its lower triangle, where \eqn{m}
    is the number of variables. If \code{FALSE},
    a vector of length \eqn{m * (m - 1) / 2} is returned, where the value for
    testing the \eqn{i}th and \eqn{j}th variable is given by the
    \eqn{j + m * (i - 1) - i * (i - 1) / 2} element of this vector.} 
}

\note{Contrary to \code{\link{chisq.test}}, currently no continuity correction is done 
  for 2 x 2 tables.
}

\value{
  If \code{compPval = FALSE}, a vector (or matrix if \code{cl} is not specified and
  \code{as.matrix = TRUE}) composed of the values of Pearson's \eqn{\chi^2}{ChiSquare}-statistic.
  Otherwise, a list consisting of
  \item{stats}{a vector (or matrix) containing the values of Pearson's \eqn{\chi^2}{ChiSquare}-statistic.}
  \item{df}{a vector (or matrix) comprising the degrees of freedom of the asymptotic
    \eqn{\chi^2}{ChiSquare}-distribution.}
  \item{rawp}{a vector (or matrix) containing the (unadjusted) p-values.} 
}


\references{
Schwender, H.\ (2007). A Note on the Simultaneous Computation of Thousands of 
Pearson's \eqn{\chi^2}{ChiSquare}-Statistics. \emph{Technical Report}, SFB 475,
Deparment of Statistics, University of Dortmund. 
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{computeContCells}}, \code{\link{computeContClass}}}

\examples{\dontrun{
# Generate an example data set consisting of 5 rows (variables)
# and 200 columns (observations) by randomly drawing integers 
# between 1 and 3.

mat <- matrix(sample(3, 1000, TRUE), 5)
rownames(mat) <- paste("SNP", 1:5, sep = "")

# For each pair of rows of mat, test if they are independent.

r1 <- rowChisqStats(mat)

# The values of Pearson's ChiSquare statistic as matrix.

r1$stats

# And the corresponding (unadjusted) p-values.

r1$rawp

# Obtain only the values of the test statistic as vector

rowChisqStats(mat, compPval = FALSE, asMatrix =FALSE)


# Generate an example data set consisting of 10 rows (variables)
# and 200 columns (observations) by randomly drawing integers 
# between 1 and 3, and a vector of class labels of length 200
# indicating that the first 100 observation belong to class 1
# and the other 100 to class 2. 

mat2 <- matrix(sample(3, 2000, TRUE), 10)
cl <- rep(1:2, e = 100)

# For each row of mat2, test if they are associated with cl.

r2 <- rowChisqStats(mat2, cl)
r2$stats

# And the results are identical to the one of chisq.test
pv <- stat <- numeric(10)
for(i in 1:10){
    tmp <- chisq.test(mat2[i,], cl)
    pv[i] <- tmp$p.value
    stat[i] <- tmp$stat
}

all.equal(r2$stats, stat)
all.equal(r2$rawp, pv)

}}

\keyword{htest}
\keyword{array}
