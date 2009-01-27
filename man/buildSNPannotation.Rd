\name{buildSNPannotation}
\alias{buildSNPannotation}

\title{Construct Annotation for Affymetrix SNP Chips}
\description{
  Constructs a data frame from a metadata package containing annotations for the SNPs
  from the corresponding Affymetrix SNP Chip.
}
\usage{
  buildSNPannotation(pkg, rs = TRUE, allele = TRUE, gene = TRUE,
    chromosome = FALSE, position = FALSE, strand = FALSE, cytoband = FALSE, 
    max.genes = 0, lib.loc = NULL, others = NULL, subset = NULL, 
    pattern = NULL, na.rm = TRUE)
}

\arguments{
  \item{pkg}{the name of the metadata package from which the data frame containing the
    annotations of the SNPs should be generated.}
  \item{rs}{should the RefSNP-ID of the SNPs be added to the data frame?}
  \item{allele}{should the two alleles of each SNP be added to the data frame?}
  \item{gene}{should the genes associated with the SNPs be added to the data frame?}
  \item{chromosome}{should the chromosome to which the respective SNP belongs be added to the data frame?}
  \item{position}{should the physical positions of the SNPs be added to the data frame?}
  \item{strand}{should the strands be added to the data frame?}
  \item{cytoband}{logical indicating whether the cytoband of each SNP is added to the data frame.}
  \item{max.genes}{integer specifying the maximum number of genes associated with the respective
    SNP that should be stored in the data frame. By default, all entries are considered. The corresponding column
    of the data frame can also be shortened afterwards using \code{\link{shortenGeneDescription}}.
    Shortened entries are marked by \code{\ldots} at the end of the entries.}
  \item{lib.loc}{the directory in which the metadata package is stored. Needs only to be specified
    if it is not stored in the usual directory of the packages.}
  \item{others}{character string or vector naming other entries of the object \code{featureSet}
    saved in \code{pkg} that should be added to the data frame.}
  \item{subset}{character string consisting of the probe set IDs of the SNPs for which the data
    frame should be generated. The data frame will contain all SNPs if \code{subset = NULL}}
  \item{pattern}{character string specifying the pattern of the probe set IDs of the SNPs for which
    the data frame should be generated. For example, \code{pattern = "SNP\%"} will lead to a
    data frame containing all SNPs whose probe set ID start with "SNP".}
  \item{na.rm}{should the rows of the data frame corresponding to SNPs specified by \code{subset}
    for which no information is available in the metadata package be removed?}
}

\value{
  A data frame composed of annotations for the SNPs for which information is available in
  the specified metadata package.
}

\author{Holger Schwender, \email{holger.schw@gmx.de}}

\seealso{\code{\link{shortenGeneDescription}}}

\keyword{misc}
