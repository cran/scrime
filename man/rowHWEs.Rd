\name{rowHWEs}
\alias{rowHWEs}

\title{Rowwise Test for Hardy-Weinberg Equilibrium}

\description{
Tests for each row of a matrix whether the Hardy-Weinberg Equilibrium holds for the SNP
represented by the row.
}
\usage{
rowHWEs(x, levels = 1:3, affy = FALSE, check = TRUE)
}

\arguments{
\item{x}{a matrix in which each row represents a SNP and each column a subject,
   where the SNPs can take the values specified by \code{levels}. NAs are allowed.}
\item{levels}{a vector of length three specifying the values with which the three
   genotypes of each SNP are represented. It is assumed that the second element of
   \code{levels} represents the heterozygous genotype, whereas the first and the
   third element represent the homozygous genotypes. Ignored if \code{affy = TRUE}.}
\item{affy}{logical specifying whether the SNPs in \code{x} are coded as in the Affymetrix
   standard output. If \code{TRUE}, \code{levels = c("AA", "AB", "BB")} will be used.}
\item{check}{should some checks be done if, e.g., other than the specified \code{levels}
   are used in \code{x}? It is highly recommended to leave \code{check = TRUE}. Setting
   \code{check = FALSE} will reduce the computation time slightly.}
}
\value{
A list containing the values of the ChiSquare statistic for testing for deviation from HWE
(\code{stats}) and the raw p-values (\code{rawp}) computed by employing the ChiSquare distribution
with 1 degree of freedom.
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\keyword{htest}
\keyword{array}

