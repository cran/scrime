\name{gknn}
\alias{gknn}

\title{Generalized k Nearest Neighbors}
\description{
  Predicts the classes of new observations with \eqn{k} Nearest Neighbors
  based on an user-specified distance measure.
}
\usage{
gknn(data, cl, newdata, nn = 5, distance = NULL, use.weights = FALSE, ...)
}

\arguments{
  \item{data}{a numeric matrix in which each row represents an observation and each column
    a variable. If \code{distance} is \code{"smc"}, \code{"cohen"} or \code{"pcc"},
    the values in \code{data} must be integers between 1 and \eqn{n_{cat}}{n.cat},
    where \eqn{n_{cat}}{n.cat} is the maximum number of levels one of the variables can
    take. Missing values are allowed.}
  \item{cl}{a numeric vector of length \code{nrow(data)} giving the class labels of
    the observations represented by the rows of \code{data}. \code{cl} must consist
    of integers between 1 and \eqn{n_{cl}}{n.cl}, where \eqn{n_{cl}}{n.cl} is the
    number of groups.}
  \item{newdata}{a numeric matrix in which each row represents a new observation for
    which the class label should be predicted and each column consists of the same
    variable as the corresponding column of \code{data}.}
  \item{nn}{an integer specifying the number of nearest neighbors used to classify
    the new observations.}
  \item{distance}{character vector naming the distance measure used to identify the
    \code{nn} nearest neighbors. Must be one of \code{"smc"}, \code{"cohen"},
    \code{"pcc"}, \code{"euclidean"}, \code{"maximum"}, \code{"manhattan"},
    \code{"canberra"}, and \code{"minkowski"}. If \code{NULL}, it is determined in
    an ad hoc way if the data seems to be categorical. If this is the case \code{distance}
    is set to \code{"smc"}. Otherwise, it is set to \code{"euclidean"}.}
  \item{use.weights}{should the votes of the nearest neighbors be weighted by the reciprocal
    of the distances to the new observation when the class of a new observation should be
    predicted?}
  \item{\dots}{further arguments for the distance measure. If, e.g., 
    \code{distance = "minkowski"}, then \code{p} can also be specified, see \code{\link{dist}}.
    If \code{distance = "pcc"}, then \code{version} can also be specified,
    see \code{\link{pcc}}.}
}

\value{
  The predicted classes of the new observations.
}

\references{
Schwender, H.\ (2007). Statistical Analysis of Genotype and Gene Expression Data.
\emph{Dissertation}, Department of Statistics, University of Dortmund.}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{knncatimpute}}, \code{\link{smc}}, \code{\link{pcc}}}

\examples{\dontrun{
# Using the example from the function knn.

library(class)
data(iris3)
train <- rbind(iris3[1:25,,1], iris3[1:25,,2], iris3[1:25,,3])
test <- rbind(iris3[26:50,,1], iris3[26:50,,2], iris3[26:50,,3])
cl <- c(rep(2, 25), rep(1, 25), rep(1, 25))

knn.out <- knn(train, test, as.factor(cl), k = 3, use.all = FALSE)
gknn.out <- gknn(train, cl, test, nn = 3)

# Both applications lead to the same predictions.

knn.out == gknn.out

# But gknn allows to use other distance measures than the Euclidean 
# distance. E.g., the Manhattan distance.

gknn(train, cl, test, nn = 3, distance = "manhattan")

}}

\keyword{classif}
