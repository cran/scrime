\name{simulateSNPs}
\alias{simulateSNPs}

\title{Simulation of SNP data}
\description{
  Simulates SNP data, where a specified proportion of cases and controls
  is explained by specified set of SNP interactions. Can also be used
  to simulate a data set with a multi-categorical response, i.e.\
  a data set in which the cases are divided into several classes (e.g.,
  different diseases or subtypes of a disease).
}
\usage{
simulateSNPs(n.obs, n.snp, vec.ia, prop.explain = 1, 
  list.ia.val = NULL, vec.ia.num = NULL, vec.cat = NULL,
  maf = c(0.1, 0.4), prob.val = rep(1/3, 3), list.equal = NULL, 
  prob.equal = 0.8, rm.redundancy = TRUE, shuffle = FALSE, 
  shuffle.obs = FALSE, rand = NA)
}

\arguments{
  \item{n.obs}{either an integer specifying the total number of 
    observations, or a vector of length 2 specifying the number 
    of cases and the number of controls. If \code{vec.cat} is specified,
    then the partitioning of the number of cases to the different
    classes can be governed by \code{vec.ia.num}. If \code{n.obs} is
    an integer, then \eqn{1 / c} of the observations will be controls and
    the remaining observations will be cases, where \eqn{c} is the total number
    of groups (including the controls).}
  \item{n.snp}{integer specifying the number of SNPs.}
  \item{vec.ia}{a vector of integers specifying the orders of the interactions
    that explain the cases. \code{c(3,1,2,3)}, e.g., means that a three-way,
    a one-way (i.e. just a SNP), a two-way, and a three-way interaction explain the cases.}
  \item{prop.explain}{either an integer or a vector of \code{length(vec.ia)} 
    specifying the proportions of cases explained by the interactions of interest
    among all observation having the interaction of interest. Must be larger than 0.5.
    E.g., \code{prop.explain = 1} means that only cases have the interactions of interest
    specified by \code{vec.ia} (and \code{list.ia.val}). E.g., \code{vec.ia = c(3, 2)}
    and \code{prop.explain = c(1, 0.8)} means that only cases have the three-way interaction of
    interest, while 80\% of the observations having the two-way interaction of interest
    are cases, and 20\% are controls.}
  \item{list.ia.val}{a list of \code{length(vec.ia)} specifying the exact interactions.
    The objects in this list must be vectors of length \code{vec.ia[i]}, and consist
    of the values 0 (for homozygous reference), 1 (heterozygous variant), or 2 (homozygous variant).
    E.g., \code{vec.ia = c(3, 2)} and \code{list.ia.val = list(c(2, 0, 1), c(0, 2))} 
    and \code{prob.equal = 1} (see also \code{list.equal}) means that
    ((SNP1 == 2) \& (SNP2 == 0) \& (SNP3 == 1)) and ((SNP4 == 0) \& (SNP5 == 2))
    are the explanatory interactions (if additionally \code{prob.equal = 1}; see also
    \code{list.equal}). If \code{NULL}, the genotypes are randomly drawn
    using the probabilities given by \code{prob.val}.}
  \item{vec.ia.num}{a vector of \code{length(vec.ia)} specifying the number of
    \emph{cases} (not observations) explained by the interactions in \code{vec.ia}.
    If \code{NULL}, all the cases are divided into \code{length(vec.ia)} groups of
    about the same size. \code{sum(vec.ia.num)} must be smaller than
    or equal to the total number of cases. Each entry of \code{vec.ia.num} must currently be >= 10.}
  \item{vec.cat}{a vector of the same length of \code{vec.ia} specifying the subclasses
    of the cases that are explained by the corresponding interaction in \code{vec.ia}. If \code{NULL},
    no subclasses will be considered. This feature is currently not fully tested. So be careful
    if specifying \code{vec.cat}.}
  \item{maf}{either an integer, or a vector of length 2 or \code{n.snp} specifying
    the minor allele frequencies. If an integer, all SNPs will have the same
    minor allele frequency. If a vector of length \code{n.snp}, each SNP will have the minor
    allele frequency specified in the corresponding entry of \code{maf}. If length 2, then
    \code{maf} is interpreted as the range of the minor allele frequencies, and for each SNP,
    a minor allele frequency will be randomly drawn from a uniform distribution with 
    the range given by \code{maf}. Note: If a SNP belongs to an explanatory interaction,
    then only the set of observations not explained by this interaction will have
    the minor allele frequency specified by \code{maf}.}
  \item{prob.val}{a vector consisting of the probabilities for drawing a 0, 1, or 2,
    if \code{list.ia.val = NULL}, i.e.\ if the genotypes of the SNPs explaining the case-control
    status should be randomly drawn. Ignored if \code{list.ia.val} is specified. By default,
    each genotype has the same probability of being drawn.}
  \item{list.equal}{list of same structure as \code{list.ia.val} containing only ones and
    zeros, where a 1 specifies the equality to the corresponding value in
    \code{list.ia.val}, and a 0 specifies the non-equality. Thus, the entries of \code{list.equal}
    specify if the corresponding SNP should be of a particular genotype (when the entry is 1)
    or should be not of this genotype (when entry is 0). If \code{NULL}, this list
    will be generated automatically using \code{prob.equal}. If, e.g., \code{vec.ia = c(3, 2)},
    \code{list.ia.val = list(c(2, 0, 1), c(0, 2))},
    and \code{list.equal = list(c(1, -1, 1), c(1, -1))}, then the explanatory interactions are
    given by ((SNP1 == 2) \& (SNP2 != 0) \& (SNP3 == 1)) and ((SNP4 == 0) \& (SNP5 != 2))}
  \item{prob.equal}{a numeric value specifying the probability that a 1 is drawn when generating
    \code{list.equal}. \code{prob.equal} is thus the probability for an equal sign.}
  \item{rm.redundancy}{should redundant SNPs be removed from the explaining interactions?
    It is possible that one specify an explaining \eqn{i}-way interaction, but an interaction
    between \eqn{(i-1)} of the variables contained in the \eqn{i}-way
    interaction already explains all the cases (and controls) that the \eqn{i}-way interaction
    should explain. In this case, the redundant SNP is removed if \code{rm.redundancy = TRUE}.}
  \item{shuffle}{logical. By default, the first \code{sum(vec.ia)} columns of the generated
    data set contain the explanatory SNPs in the same order as they appear in this data set.
    If \code{TRUE}, this order will be shuffled.}
  \item{shuffle.obs}{should the observations be shuffled?}
  \item{rand}{integer. Sets the random number generator in a reproducible state.}
}

\note{Currently, the genotypes of all SNPs are simulated independently from each other
  (except for the SNPs that belong to the same explanatory interaction).}
  
\value{
  An object of class \code{simulatedSNPs} composed of
  \item{data}{a matrix with \code{n.obs} rows and \code{n.snp} columns containing
    the SNP data.}
  \item{cl}{a vector of length \code{n.obs} comprising the case-control status of the
    observations.}
  \item{tab.explain}{a table naming the explanatory interactions and the numbers
    of cases and controls explained by them.}
  \item{ia}{character vector naming the interactions.}
  \item{maf}{vector of length \code{n.snp} containing the minor allele frequencies.}
}

\author{Holger Schwender \email{holger.schwender@udo.edu}}

\seealso{\code{\link{simulateSNPglm}}, \code{\link{simulateSNPcatResponse}}}

\examples{\dontrun{
# Simulate a data set containing 2000 observations (1000 cases
# and 1000 controls) and 50 SNPs, where one three-way and two 
# two-way interactions are chosen randomly to be explanatory 
# for the case-control status.

sim1 <- simulateSNPs(2000, 50, c(3, 2, 2))
sim1

# Simulate data of 1200 cases and 800 controls for 50 SNPs, 
# where 90% of the observations showing a randomly chosen 
# three-way interaction are cases, and 95% of the observations 
# showing a randomly chosen two-way interactions are cases.

sim2 <- simulateSNPs(c(1200, 800), 50, c(3, 2), 
   prop.explain = c(0.9, 0.95))
sim2

# Simulate a data set consisting of 1000 observations and 50 SNPs,
# where the minor allele frequency of each SNP is 0.25, and
# the interactions 
# ((SNP1 == 2) & (SNP2 != 0) & (SNP3 == 1))   and 
# ((SNP4 == 0) & (SNP5 != 2))
# are explanatory for 200 and 250 of the 500 cases, respectively,
# and for none of the 500 controls.

list1 <- list(c(2, 0, 1), c(0, 2))
list2 <- list(c(1, 0, 1), c(1, 0))
sim3 <- simulateSNPs(1000, 50, c(3, 2), list.ia.val = list1,
    list.equal = list2, vec.ia.num = c(200, 250), maf = 0.25)

}}

\keyword{datagen}
