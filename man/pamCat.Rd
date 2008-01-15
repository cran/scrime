\name{pamCat}
\alias{pamCat}
\alias{print.pamCat}

\title{Prediction Analysis of Categorical Data}
\description{
 Performs a Prediction Analysis of Categorical Data.
}
\usage{
pamCat(data, cl, theta = NULL, n.theta = 10, newdata = NULL, newcl = NULL)
}

\arguments{
  \item{data}{a numeric matrix composed of the integers between 1 and \eqn{n_{cat}}{n.cat},
    where \eqn{n_{cat}}{n.cat} is the number of levels each of the variables represented
    by the rows of \code{data} must take. No missing values allowed.}
  \item{cl}{a numeric vector of length \code{ncol(data)} comprising the class labels of
    the observations represented by the columns of \code{data}. \code{cl} must consist
    of the integers between 1 and \eqn{n_{cl}}{n.cl}, where \eqn{n_{cl}}{n.cl} is the number
    of classes.}
  \item{theta}{a numeric vector consisting of the strictly positive values of the shrinkage parameter used
    in the Prediction Analysis. If \code{NULL}, a vector consisting of \code{n.theta} values for
    the shrinkage parameter are determined automatically.}
  \item{n.theta}{an integer specifying the number of values for the shrinkage parameter of the
    Prediction Analysis. Ignored if \code{theta} is specified.}
  \item{newdata}{a numeric matrix composed of the integers between 1 and \eqn{n_{cat}}{n.cat}.
    Must have the same number of rows as \code{data}, and each row of \code{newdata} must contain
    the same variable as the corresponding row of \code{data}. \code{newdata} is employed to
    compute the misclassification rates of the Prediction Analysis for the given values of the
    shrinkage parameter. If \code{NULL}, \code{data} is used to determine the misclassification rates.}
  \item{newcl}{a numeric vector of length \code{ncol(newdata)} that consists of integers between
    1 and \eqn{n_{cl}}{n.cl}, and specifies the class labels of the observations in \code{newdata}.
    Must be specified, if \code{newdata} is specified.}
}

\value{
  An object of class \code{pamCat} composed of
  \item{mat.chisq}{a matrix with \eqn{m} rows and \eqn{n_{cl}}{n.cl} columns consisting of the
    classwise values of Pearson's \eqn{\chi^2}{ChiSquare} statistic for each of the \eqn{m} variables.}
  \item{mat.obs}{a matrix with \eqn{m} rows and \eqn{n_{cat} * n_{cl}}{n.cat * n.cl} columns
    in which each row shows a contingency table between the corresponding variable and \code{cl}.}
  \item{mat.exp}{a matrix of the same size as \code{mat.obs} containing the numbers of observations
    expected under the null hypothesis of an association between the respective variable and \code{cl}.}
  \item{mat.theta}{a data frame consisting of the numbers of variables used in the classification
    of the observations in \code{newdata} and the corresponding misclassification rates for a set of values of
    the shrinkage parameter \eqn{\theta}{theta}.}
  \item{tab.cl}{a table summarizing the values of the response, i.e.\ the class labels.}
  \item{n.cat}{\eqn{n_{cat}}{n.cat}.}
}

\references{
Schwender, H.\ (2007). Statistical Analysis of Genotype and Gene Expression Data.
\emph{Dissertation}, Department of Statistics, University of Dortmund.
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{predict.pamCat}}}

\examples{\dontrun{
# Generate a data set consisting of 2000 rows (variables) and 50 columns.
# Assume that the first 25 observations belong to class 1, and the other
# 50 observations to class 2.

mat <- matrix(sample(3, 100000, TRUE), 2000)
rownames(mat) <- paste("SNP", 1:2000, sep = "")
cl <- rep(1:2, e = 25)

# Apply PAM for categorical data to this matrix, and compute the
# misclassification rate on the training set, i.e. on mat.

pam.out <- pamCat(mat, cl)
pam.out

# Now generate a new data set consisting of 20 observations, 
# and predict the classes of these observations using the
# value of theta that has led to the smallest misclassification
# rate in pam.out.

mat2 <- matrix(sample(3, 40000, TRUE), 2000)
rownames(mat2) <- paste("SNP", 1:2000, sep = "")
predict(pam.out, mat2)

# Let's assume that the predicted classes are the real classes
# of the observations. Then, mat2 can also be used in pamCat
# to compute the misclassification rate. 

cl2 <- predict(pam.out, mat2)
pamCat(mat, cl, newdata = mat2, newcl = cl2)

}}


\keyword{classif}
