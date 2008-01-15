\name{recodeAffySNP}
\alias{recodeAffySNP}

\title{Recoding of Affymetrix SNP Values}
\description{
  Recodes the values used on Affymetrix SNP chips to code the genotypes to
  other values -- required, e.g., by other functions of this package.
}
\usage{
recodeAffySNP(mat, refAA = FALSE, geno = 1:3)
}

\arguments{
  \item{mat}{a matrix or data frame consisting of the character strings \code{"AA"}, \code{"AB"},
    \code{"BB"}. Missing values can be coded by either \code{"NN"} or \code{NA}. Each
    row is assumed to correspond to a variable, and each column to a microarray.}
  \item{refAA}{codes \code{"AA"} always for the homozygous reference gentoype? If
   \code{TRUE}, \code{"AA"} is always replaced by \code{geno[1]}, and \code{"BB"} by
   \code{geno[3]}. If \code{FALSE}, it is evaluated rowwise whether \code{"AA"} or
   \code{"BB"} occurs more often, and the more frequently occuring value is set to \code{geno[1]}.}
  \item{geno}{a numeric or character vector of length 3 giving the three values that
    should be used to recode the genotypes. By default, \code{geno = 1:3} which is the
    coding, e.g., required by \code{\link{rowChisqStats}} or \code{\link{pamCat}}.}
}
\value{
  A matrix of the same size as \code{mat} containing the recoded genotypes. (Missing values are
  coded by \code{NA}.)
}

\seealso{\code{\link{recodeSNPs}}, \code{\link{snp2bin}}}

\examples{\dontrun{
# Generate a sample data set consisting of 10 rows and 12 columns,
# and randomly replace 5 of the values by "NN".

mat <- matrix("", 10, 12)
mat[1:5,] <- sample(c("AA", "AB", "BB"), 60, TRUE,
   prob = c(0.49, 0.42, 0.09))
mat[6:10,] <- sample(c("AA", "AB", "BB"), 60, TRUE,
   prob = c(0.09, 0.42, 0.49))
mat[sample(120, 5)] <- "NN"
mat

# Recode the SNPs.

recodeAffySNP(mat)

# Recode the SNPs assuming that "A" is always the major allele.

recodeAffySNP(mat, refAA = TRUE)

}}




\keyword{manip}
\keyword{array}
