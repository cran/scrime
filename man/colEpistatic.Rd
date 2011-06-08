\name{colEpistatic}
\alias{colEpistatic}
\alias{rowEpistatic}
\alias{print.colEpi}
\title{Cordell's Test for Epistatic Interactions
}
\description{Performs the likelihood ratio test for epistatic interactions proposed by Cordell (2002) for case-control data, where
  \code{colEpistatic} assumes that each column represents a SNP, and \code{rowEpistatic} assumes that each row represents a SNP.
}
\usage{
colEpistatic(mat.snp, cl, genes = NULL, warnError = TRUE)
rowEpistatic(mat.snp, cl, genes = NULL, warnError = TRUE)

\method{print}{colEpi}(x, top = 5, digits = 4, ...)
}

\arguments{
  \item{mat.snp}{a matrix containing genotype data, where the genotypes of each SNPs need to be coded by the number of minor alleles,
    i.e. 0, 1, and 2. Missing values are allowed. For \code{colEpistatic}, each column of \code{mat.snp} is assumed to represent a SNP, and
    each row a subject. For \code{rowEpistatic}, each row of \code{mat.snp} is assumed to represent a subject, and each column a SNP. 
}
  \item{cl}{a numeric vector of ones and zeros specifying which of the subjects in \code{mat.snp} are cases (coded by 1) and which are
    controls (coded by 0).
}
  \item{genes}{a character vector containing the names of the genes (or, e.g., LD-blocks or pathways) to which the SNPs belong. If specified, 
    only the two-way interactions between SNPs from different genetic sets (e.g., genes, LD-blocks, or pathways) are tested. 
    If \code{NULL}, all two-way interactions between all possible pairs of SNPs are tested.
}
  \item{warnError}{logical indicating whether the statistics for the gTDT for pairs of SNPs should be returned as \code{NA} if the fitting 
    of one of the logistic regression models with \code{\link[stats]{glm}} returns warning. If \code{warnError = FALSE}, 
    these warnings will be ignored, and the gTDT will be performed also for these SNP pairs.
}
  \item{x}{an object of class \code{colEpi}, i.e. the output of either \code{colEpistatic} or \code{rowEpistatic}.
}
 \item{top}{number of interactions that should be printed. If \code{top} is set to \code{0} or the number of interactions is smaller than 
    or equal to \code{top}, then the statistics for all interactions are printed in the order of their computation. Otherwise, 
    the \code{top} interactions with the smallest p-values are printed (sorted by their p-values).
}
  \item{digits}{number of digits that should be printed.
}
  \item{...}{ignored.
}
}

\value{
An object of class \code{colEpi} consisting of
\item{ll.main}{a numeric vector containing the values of the maximized loglikelihoods of the logistic regression models 
   considering only main effects,}
\item{ll.full}{a numeric vector containing the values the maximized loglikelihoods of the logistic regression models additionally 
   containing interaction terms,}
\item{stat}{a vector comprising the values of the test statistic,}
\item{pval}{a vector comprising the corresponding p-values,}
\item{genes}{if \code{genes} has been specified, a vector composed of the gene pairs to which the SNP pairs belong,}
\item{vec.error}{if \code{warnError = TRUE}, a list in which an entry corresponding to a SNP pair for which the fitting of (at least)
  one of the logistic regression models led to a warning shows this warning message of \code{\link[stats]{glm}}, and all other entries
  are \code{NULL}.}
}

\references{
Cordell, H. J. (2002). Epistasis: What it Means, what it Doesn't mean, and Statistical Methods to Detect it in Humans.
Human Molecular Genetics, 11, 2463-2468.
}

\author{
Holger Schwender, \email{holger.schwender@udo.edu}
}

\seealso{
\code{\link[trio]{colTDT2way}}
}

\keyword{htest}
\keyword{models}
\keyword{array}
