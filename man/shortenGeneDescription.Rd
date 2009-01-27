\name{shortenGeneDescription}
\alias{shortenGeneDescription}

\title{Shorten the Gene Description}
\description{
  Shortens the entries of the column of a data frame containing the genes
  associated with the SNPs for which the data frame comprises annotations.
  Typically used in combination with, i.e.\ either within or after an application of,
  \code{\link{buildSNPannotation}.} 
}
\usage{
  shortenGeneDescription(dat, colname = "Gene", max.length = 2, 
     sep = "///", add.ldots = TRUE)
}

\arguments{
  \item{dat}{a data frame. Typically, the output of \code{\link{buildSNPannotation}}.}
  \item{colname}{character string comprising the name of the column of \code{dat}
    containing the gene description.}
  \item{max.length}{integer specifying the maximum number of genes associated with the respective
    SNP that should be stored in the data frame. By default, the first two genes are retained. 
    Shortened entries are marked by \code{\ldots} at the end of the entries, when \code{add.ldots = TRUE}.}
  \item{sep}{character string specifying the separation symbol between the different genes.}
  \item{add.ldots}{should \code{\ldots} be added at the entries which are shortened?}
}

\value{
 The same data frame as \code{dat} with shortened entries in the column of \code{dat} named \code{colname}.
}
\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{\code{\link{buildSNPannotation}}}
   
\keyword{utilities}
