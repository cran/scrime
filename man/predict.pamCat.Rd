\name{predict.pamCat}
\alias{predict.pamCat}

\title{Predict Method for pamCat Objects}
\description{
  Predicts the classes of new observations based on a Prediction Analysis of Categorical Data.
}
\usage{
\method{predict}{pamCat}(object, newdata, theta = NULL, add.nvar = FALSE,
   type = c("class", "prob"), ...)
}

\arguments{
  \item{object}{an object of class \code{pamCat}.}
  \item{newdata}{a numeric matrix consisting of the integers between 1 and \eqn{n_{cat}}{n.cat}, where
    \eqn{n_{cat}}{n.cat} is the number of levels each of the variables in \code{newdata} must take. 
    Each row of \code{newdata} must represent the same variable as the corresponding row in
    the matrix \code{data} used to produce \code{object}.
    Each column corresponds to an observation for which the class should be predicted.}
  \item{theta}{a strictly positive numeric value specifying the value of the shrinkage parameter
    of the Prediction Analysis that should be used in the class prediction.
    If \code{NULL}, then the value of \code{theta} will be used that
    has led to the smallest misclassification rate in the initial Prediction Analysis generating \code{object}.}
  \item{add.nvar}{should the number of variables used in the class prediction be added to the output?}
  \item{type}{either \code{"class"} or \code{"prob"}. If \code{"class"}, then the predicted classes will be returned.
    Otherwise, the probabilities for the classes are returned.}
  \item{\dots}{Ignored.}
}
\value{
  If \code{add.nvar = FALSE}, the predicted classes or the class probabilities (depending on \code{type}). 
  Otherwise, a list consisting of
  \item{pred}{a vector or matrix containing the predicted classes or the class probabilities, respectively.}
  \item{n.var}{the number of variables used in the prediction.}
}

\references{
Schwender, H.\ (2007). Statistical Analysis of Genotype and Gene Expression Data.
\emph{Dissertation}, Department of Statistics, University of Dortmund.
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{pamCat}}}

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

# Another theta, say theta = 4, can also be specified.

predict(pam.out, mat2, theta = 4)

# The class probabilities for each observation can be obtained by

predict(pam.out, mat2, theta = 4, type = "prob") 

}}

\keyword{classif}
