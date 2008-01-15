\name{snp2bin}
\alias{snp2bin}

\title{Transformation of SNPs to Binary Variables}
\description{
  Transforms SNPs to binary variables.
}
\usage{
snp2bin(mat, domrec = TRUE, refAA = FALSE, snp.in.col = TRUE, 
   monomorph = 0)
}

\arguments{
  \item{mat}{a matrix or data frame in which the genotypes of all SNPs are coded
    either by \code{0}, \code{1} and \code{2},
    or by \code{1}, \code{2} and \code{3}, or by \code{"AA"}, \code{"AB"} and \code{"BB"}.
    Missing values are allowed. In the latter coding not only \code{NA}, but also \code{"NN"}
    is allowed for specifying missing values. Using the former two codings it is assumed
    that the smallest value codes the homozygous reference genotype, the second value the
    heterozygous genotype, and the largest value the homozygous variant genotype. For the third
    coding, see \code{refAA}.}
  \item{domrec}{should each SNP be coded by two dummy variables from which one codes for
    a recessive, and the other for a dominant effect? If \code{TRUE}, then the first
    binary variable is set to 1 if the SNP is of the heterozygous or the homozygous variant
    genotype, and the second dummy variable is set to 1 if the SNP is of the homozygous
    variant genotype. If \code{FALSE}, three dummy variables are used and each of the three
    genotypes of a SNP is coded by one of these binary variables.}
  \item{refAA}{codes \code{"AA"} always for the homozygous reference gentoype? Only considered
   if the SNPs are coded by \code{"AA"}, \code{"AB"} and \code{"BB"}.
   If \code{FALSE}, it is evaluated SNPwise whether \code{"AA"} or
   \code{"BB"} occurs more often, and the more frequently occuring value is assumed to be
   the homozygous reference genotype.}
  \item{snp.in.col}{does each column of \code{mat} correspond to a SNP (and each row to an observation)?
    If \code{FALSE}, it is assumed that each row represents a SNP, and each column an observation.}
  \item{monomorph}{a non-negative number. If a dummy variable contains \code{monomorph} or less
    values that differ from the more frequent value of this variable, then the variable is
    removed from the data set.}
}
\value{
   A matrix containing the binary dummy variables.
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{recodeSNPs}}, \code{\link{recodeAffySNP}}}

\examples{\dontrun{
# Generate an example data set consisting of 10 rows (observations)
# and 5 columns (SNPs).

mat <- matrix(sample(3, 50, TRUE), 10)
colnames(mat) <- paste("SNP", 1:5, sep = "")

# Transform each SNP into two dummy variables, one that codes for
# a recessive effect and one that codes for a dominant effect.

snp2bin(mat)

# Transform each SNP into three dummy variables, where each of
# these variables codes for one of the three genotypes.

snp2bin(mat, domrec = FALSE)
}} 

\keyword{array}
\keyword{manip}
