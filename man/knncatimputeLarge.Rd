\name{knncatimputeLarge}
\alias{knncatimputeLarge}

\title{Missing Value Imputation with kNN for High-Dimensional Data}
\description{
  Imputes missing values in a high-dimensional matrix composed of categorical variables
  using \eqn{k} Nearest Neighbors.
}
\usage{
knncatimputeLarge(data, mat.na = NULL, fac = NULL, fac.na = NULL,
   nn = 3, distance = c("smc", "cohen", "snp1norm", "pcc"), 
   n.num = 100, use.weights = TRUE, verbose = FALSE)
}

\arguments{
  \item{data}{a numeric matrix consisting of integers between 1 and \eqn{n_{cat}}{n.cat},
    where \eqn{n_{cat}}{n.cat} is maximum number of levels the categorical variables
    can take. If \code{mat.na} is specified, \code{data} is assumed to contain only
    non-missing data, and the rows of \code{data} are used to impute the missing values
    in \code{mat.na}. Otherwise, \code{data} is also allowed to contain missing values,
    and the missing values in the rows of \code{data} are imputed by employing the rows of
    \code{data} showing no missing values.
    
    Each row of \code{data} represents one of the objects that should be used to identify the
    \eqn{k} nearest neighbors, i.e.\ if the \eqn{k} nearest variables should be used to
    replace the missing values, then each row must represent one of the variables. If the
    \eqn{k} nearest observations should be used to impute the missing values, then each
    row must correspond to one of the observations.}
  \item{mat.na}{a numeric matrix containing missing values. Must have the same number of
    columns as \code{data}. All non-missing values must be integers between 1 and
    \eqn{n_{cat}}{n.cat}. If \code{NULL}, \code{data} is assumed to also contain the
    rows with missing values.}
  \item{fac}{a numeric or character vector of length \code{nrow{data}} specifying the values of
    a factor used to split \code{data} into subsets. If, e.g., the values of \code{fac}
    are given by the chromosomes to which the SNPs represented by the rows of \code{data}
    belong, then \eqn{k} nearest neighbors is applied chromosomewise to the missing values
    in \code{mat.na} (or \code{data}). If \code{NULL}, no such splitting is done. Must be
    specified, if \code{fac.na} is specified.}
  \item{fac.na}{a numeric or character vector of length \code{nrow{mat.na}} specifying the values
    of a factor by which \code{mat.na} is split into subsets. Each possible value of \code{fac.na}
    must be at least \code{nn} times in \code{fac}. Must be specified, if \code{fac} and \code{mat.na}
    is specified. If both \code{fac} and \code{fac.na} are \code{NULL}, then no splitting is done.}
  \item{nn}{an integer specifying \eqn{k}, i.e.\ the number of nearest neighbors, used to impute
    the missing values.}
  \item{distance}{character string naming the distance measure used in \eqn{k} Nearest Neighbors.
    Must be either \code{"smc"} (default), \code{"cohen"}, \code{"snp1norm"} (which denotes the
    Manhattan distance for SNPs), or \code{"pcc"}.}
  \item{n.num}{an integer giving the number of rows of \code{mat.na} considered simultaneously
    when replacing the missing values in \code{mat.na}.}
  \item{use.weights}{should weighted \eqn{k} nearest neighbors be used to impute the missing values?
    If \code{TRUE}, the votes of the nearest neighbors are weighted by the reciprocal of their
    distances to the variable (or observation) whose missing values are imputed.}
  \item{verbose}{should more information about the progress of the imputation be printed?}
}

\note{
  While in \code{\link{knncatimpute}} all variable/rows are considered when replacing
  missing values, \code{knncatimputeLarge} only considers the rows with no missing values
  when searching for the \eqn{k} nearest neighbors.
}

\value{
  If \code{mat.na = NULL}, then a matrix of the same size as \code{data} in which the missing
  values have been replaced. If \code{mat.na} has been specified, then a matrix of the same size as
  \code{mat.na} in which the missing values have been replaced.
}

\references{
Schwender, H. and Ickstadt, K.\ (2008). Imputing Missing Genotypes with \eqn{k} Nearest Neighbors.
\emph{Technical Report}, SFB 475, Department of Statistics, University of Dortmund. Appears soon.
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{knncatimpute}}, \code{\link{gknn}}, \code{\link{smc}}, \code{\link{pcc}}}

\examples{\dontrun{
# Generate a data set consisting of 100 columns and 2000 rows (actually,
# knncatimputeLarge is made for much larger data sets), where the values
# are randomly drawn from the integers 1, 2, and 3.
# Afterwards, remove 200 of the observations randomly.

mat <- matrix(sample(3, 200000, TRUE), 2000)
mat[sample(200000, 20)] <- NA

# Apply knncatimputeLarge to mat to remove the missing values.

mat2 <- knncatimputeLarge(mat)
sum(is.na(mat))
sum(is.na(mat2))

# Now assume that the first 100 rows belong to SNPs from chromosome 1,
# the second 100 rows to SNPs from chromosome 2, and so on.

chromosome <- rep(1:20, e = 100)

# Apply knncatimputeLarge to mat chromosomewise, i.e. only consider
# the SNPs that belong to the same chromosome when replacing missing
# genotypes.

mat4 <- knncatimputeLarge(mat, fac = chromosome)

}}


\keyword{classif}
\keyword{NA}
