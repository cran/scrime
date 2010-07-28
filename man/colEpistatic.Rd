\name{colEpistatic}
\alias{colEpistatic}
\alias{print.colEpi}
\title{Testing Epistatic Interactions
}
\description{
Tests epistatic two-way interactions in case-control study based on the lielihood ratio test proposed by Cordell (2002).
}
\usage{
colEpistatic(X, y)

\method{print}{colEpi}(x, top = 5, digits = 4, ...)
}
\arguments{
  \item{X}{a matrix containing the genotypes coded by 0, 1, and 2, i.e. by the number of minor alleles.
    Each pair of columns of \code{X} is tested with the likelihood ratio test of Cordell (2002). Thus, each
    column is assumed to represent one of the SNPs. This is in contrast to most other functions in this package
    in which the rows usually represent the SNPs.
}
  \item{y}{a vector of 0's (coding for controls) and 1's (coding for cases)
}
  \item{x}{the output of \code{colEpistatic}
}
  \item{top}{number of interactions that should be printed. If the number of interactions is smaller than or equal to 
    \code{top}, then the statistics for all interactions are printed in the order of their computation. Otherwise, 
    the \code{top} interactions with the smallest p-values are printed.
}
  \item{digits}{number of digits that should be printed.}
  \item{...}{ignored.}
}
\value{
An object of class \code{colEpi} consisting of

\item{ll.main}{the maximized loglikelihoods of the logistic regression models considering only main effects,}
\item{ll.full}{the maximized loglikelihoods of the logistic regression models additionally containing interaction terms,}
\item{stat}{the values of the test statistic,}
\item{pval}{the corresponding p-values.}
}

\references{
Cordell, H. J. (2002). Epistasis: What it Means, what it Doesn't mean, and Statistical Methods to Detect it in Humans.
Human Molecular Genetics, 11, 2463-2468.
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\keyword{htest}
\keyword{models}
\keyword{array}
