\name{rowMAFs}
\alias{rowMAFs}

\title{
Rowwise Minor Allele Frequency
}
\description{
Computes for each SNP represented by a row of a matrix the frequency of the minor allele.
}
\usage{
rowMAFs(x, check = TRUE)
}

\arguments{
\item{x}{a matrix in which each row represents a SNP and each column a subject,
   where the genotypes of each SNP are coded by 1 (for the homozygous reference genotype),
   2 (heterozygous), and 3 (homozygous variant). NAs are also allowed.}
\item{check}{should it be checked if the matrix contains values differing from 1, 2, and 3?
   It is highly recommended to leave \code{check = TRUE}. Setting \code{check = FALSE}
   reduces the computation time only slightly.}
}
\value{a vector containing the minor allele frequency of the SNPs represented by \code{x}.
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\keyword{array}
\keyword{manip}
