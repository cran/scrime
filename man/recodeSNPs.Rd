\name{recodeSNPs}
\alias{recodeSNPs}

\title{Recoding of SNP Values}
\description{
  Recodes the values used to specify the genotypes of the SNPs to
  other values. Such a recoding might be required to use other
  functions contained in this package.
}
\usage{
recodeSNPs(mat, first.ref = FALSE, geno = 1:3, snp.in.col = FALSE)
}

\arguments{
  \item{mat}{a matrix or data frame consisting of character strings of length 2 that
    specify the genotypes of the SNPs. Each of these character strings
    must be a combination of the letters A, T, C, and G. Missing values can
    be specified by \code{"NN"} or \code{NA}. Depending on
    \code{snp.in.col} it is assumed that each row of \code{mat} represents
    a SNP and each column a variable (\code{snp.in.col = FALSE}), or vice versa.}
  \item{first.ref}{does the first letter in the string coding the heterozygous
    genotype always stands for the more frequent allele? E.g., codes \code{"CC"}
    for the homozygous reference genotype if the genotypes
    of a SNP are coded by \code{"CC"}, \code{"CG"} and \code{"GG"}? If \code{TRUE},
    the value made up only of this first letter is set to \code{geno[1]}, and the
    value made up only of the second letter is set to \code{geno[3]}. If \code{FALSE},
    it is evaluated rowwise which of the homozygous genotypes has the higher frequency
    and the more often occuring value is set to \code{geno[1]}, and the other to \code{geno[3]}.}
  \item{geno}{a numeric or character vector of length 3 giving the three values that
    should be used to recode the genotypes. By default, \code{geno = 1:3} which is the
    coding, e.g., required by \code{\link{rowChisqStats}} or \code{\link{pamCat}}.}
  \item{snp.in.col}{does each column of \code{mat} correspond to a SNP (and each row to an array)?
    If \code{FALSE}, it is assumed that each row represents a SNP, and each column an array.}
}
\value{
  A matrix of the same size as \code{mat} containing the recoded genotypes. (Missing values are
  coded by \code{NA}).
}

\seealso{\code{\link{recodeAffySNP}}, \code{\link{snp2bin}}}

\examples{\dontrun{
# Generate an example data set consisting of 5 rows and 12 columns,
# where it is assumed that each row corresponds to a SNP.

mat <- matrix("", 10, 12)
mat[c(1, 4, 6),] <- sample(c("AA", "AT", "TT"), 18, TRUE)
mat[c(2, 3, 10),] <- sample(c("CC", "CG", "GG"), 18, TRUE)
mat[c(5, 8),] <- sample(c("GG", "GT", "TT"), 12, TRUE)
mat[c(7, 9),] <- sample(c("AA", "AC", "CC"), 12, TRUE)
mat

# Recode the SNPs

recodeSNPs(mat)

# Recode the SNPs by assuming that the first letter in
# the heterogyzous genotype refers to the major allele.

recodeSNPs(mat, first.ref = TRUE)

}}

\keyword{manip}
\keyword{array}
