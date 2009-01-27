\name{rowFreqs}
\alias{rowFreqs}

\title{Rowwise Frequencies}
\description{
  Computes the frequencies of the levels that the categorical variables
  in a matrix show.
}
\usage{
rowFreqs(x, levels = 1:3, divide.by.n = FALSE, affy = FALSE,
   includeNA = FALSE, useNN = c("not", "only", "also"), check = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a matrix in which each row represents a categorical variable (e.g., a SNP)
     and each column an observation, where the variables are assumed to show the
     levels specified by \code{levels}. Missing values are allowed in \code{x}.}
  \item{levels}{vector specifying the levels that the categorical variables in \code{x}
     show. Ignored if \code{affy = TRUE}.}
  \item{divide.by.n}{should the numbers of observations showing the respective levels
     be divided by the total number of observations, i.e.\ by \code{ncol(x)}? If \code{FALSE},
     these numbers are divided by the number of non-missing values of the respective variable.
     Ignored if \code{includeNA = TRUE}.}
  \item{affy}{logical specifying whether the SNPs in \code{x} are coded in the Affymetrix
     standard way. If \code{TRUE}, \code{levels = c("AA", "AB", "BB")} and \code{useNN = "also"}
     will be used (the latter only when \code{includeNA = TRUE}).}
  \item{includeNA}{should a column be added to the output matrix containing the number of
     missing values for each variable?}
  \item{useNN}{character specifying whether missing values can also be coded by \code{"NN"}.
     If \code{useNN = "not"} (default), missing values are assumed to be coded only by \code{NA}.
     If \code{"only"}, then missing values are assumed to be coded only by \code{"NN"} (and not
     by \code{NA}. If \code{"both"}, both \code{"NN"} and \code{NA} are considered. Ignored
     if \code{affy = TRUE}.}
  \item{check}{should it be checked whether some of the variables show other levels than the one
    specified by \code{levels}?}
}

\value{
  A matrix with the same number of rows as \code{x} containing for each variable the numbers
  of observations showing the levels specified by \code{levels}.
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{rowTables}}}

\examples{\dontrun{
# Generate a matrix containing data for 10 categorical 
# variables with levels 1, 2, 3.

mat <- matrix(sample(3, 500, TRUE), 10)

rowFreqs(mat)

# leads to the same results as 

rowTables(mat) / ncol(mat)

# If mat contains missing values

mat[sample(500, 20)] <- NA

# then

rowFreqs(mat)

# leads to the same result as

rowTables(mat) / rowSums(!is.na(mat))

}}

\keyword{array}
\keyword{manip}
