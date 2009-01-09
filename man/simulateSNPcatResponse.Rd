\name{simulateSNPcatResponse}
\alias{simulateSNPcatResponse}
\alias{print.simSNPcatResponse}

\title{Simulation of SNP Data with Categorical Response}
\description{
  Simulates SNP data. Interactions of some of the simulated SNPs
  are then used to specify a categorical response
  by level-wise logistic regression models.
}
\usage{
simulateSNPcatResponse(n.obs = 1000, n.snp = 50, list.ia = NULL,
   list.snp = NULL, beta0 = -0.5, beta = 1.5, maf = 0.25,
   sample.y = TRUE, rand = NA)
  
\method{print}{simSNPcatResponse}(x, justify = c("left", "right"), spaces = 2, ...)
}

\arguments{
  \item{n.obs}{number of observations that should be generated.}
  \item{n.snp}{number of SNPs that should be generated.}
  \item{list.ia}{a list consisting of \eqn{n_{cat}}{n.cat} objects, where
     \eqn{n_{cat}}{n.cat} is the number of levels the response should have. If one
     interaction of SNPs should be explanatory for a specific level of the response,
     then the corresponding object in \code{list.ia} must be a numeric vector specifying
     the genotypes of the interacting SNPs by the integers -3, -2, -1, 1, 2, or 3, where
     1 codes for the homozygous reference genotype, 2 for the heterozygous genotype, and
     3 for the homozygous variant genotype, and a minus before these numbers means that
     the corresponding SNP should be not of this genotype. If more than one interaction
     should be explanatory for a specific category, then the corresponding object of
     \code{list.ia} must be a list containing one numeric vector composed of the integers
     -3, -2, -1, 1, 2, and 3 for each of the interactions.
    
     If, e.g., one of the
     vectors is given by \code{c(1, -1, -3)} and the corresponding vector in
     \code{list.snp} is \code{c(5, 7, 8)}, then the corresponding interaction
     explanatory for a level of the response is given by
      
     \code{(SNP5 == 1) & (SNP7 != 1) & (SNP8 != 3)}.
     
     For more details, see Details. Must be specified if \code{list.snp} is specified.
     If both \code{list.ia} and \code{list.snp} are \code{NULL}, then the interactions
     shown in the Details section are used. }
  \item{list.snp}{a list consisting of numeric vectors (if one interaction should be explanatory
     for a level of the response) or lists of numeric vectors (if there should be more than one
     explanatory interaction) specifying the SNPs that compose
     the interactions. \code{list.snp} must have the same structure as \code{list.ia}, and each
     entry of \code{list.snp} must be an integer between 1 and \code{n.snp}.
     If \code{list.ia} is specified but
     not \code{list.snp}, then the first \eqn{n} SNPs are used to generate the interactions,
     where \eqn{n} is the total number of values in \code{list.ia}. For the case that both
     \code{list.ia} and \code{list.snp} are not specified, see Details.}
  \item{beta0}{a numeric value or vector of \code{length(list.ia)} specifying the intercept of
     the logistic regression models.}
  \item{beta}{either a non-negative numeric value or a list of non-negative numeric values specifying
    the parameters in the logistic regression models. If a numeric value, all parameters (except for
    the intercept) in all logistic regression models will be equal to this value. If a list, then
    this list must have the same length as \code{list.ia}, and each object must consist of as many
    numeric values as interactions are specified by the corresponding object in \code{list.ia}.}
  \item{maf}{either an integer, or a vector of length 2 or \code{n.snp} specifying the
    minor allele frequency. If an integer, all the SNPs will have the same
    minor allele frequency. If a vector of length \code{n.snp}, each SNP will have the minor
    allele frequency specified in the corresponding entry of \code{maf}. If length 2, then
    \code{maf} is interpreted as the range of the minor allele frequencies, and for each SNP,
    a minor allele frequency will be randomly drawn from a uniform distribution with
    the range given by \code{maf}.}
  \item{sample.y}{should the values of the response be randomly drawn using the probabilities
    determined by the logistic regression models? If \code{FALSE}, then for each of the \code{n.obs}
    observations, the value of the response is given by the level exhibiting the largest probability
    at this observation.}
  \item{rand}{a numeric value for setting the random number generator in a reproducible state.}
  \item{x}{the output of \code{simulateSNPcatResponse}}
  \item{justify}{a character string specifying whether the column of the summarizing table that
    names the explanatory interactions should be \code{"left"}- or \code{"right"}-adjusted.}
  \item{spaces}{integer specifying the distance from the left end of the column mentioned in \code{justify}
    to the position at which the column name is presented.}
  \item{\dots}{ignored.}
}
\details{
  \code{simulateSNPcatResponse} first simulates a matrix consisting of \code{n.obs} observations and \code{n.snp}
  SNPs, where the minor allele frequencies of these SNPs are given by \code{maf}.
  
  Note that all SNPs are currently simulated independently of each other such that they are unlinked.
  Moreover, an observation is currently not allowed to have genotypes/interactions that are explanatory
  for more than one of the levels of the response. If, e.g., the response has three categories, then
  an observation can either exhibit one (or more) of the genotypes explaining the first level, or
  one (or more) of the genotypes explanatory for the second level, or one (or more) of the genotypes
  explaining the third level, or none of these genotypes.

  Afterwards, the response is determined using level-wise logistic regression models specified by \code{list.ia},
  \code{list.snp}, \code{beta0} and \code{beta}.
   
  By default, i.e.\ if both \code{list.ia} and \code{list.snp} are \code{NULL}, \code{list.ia} is set
  to
  
  \code{list(c(-1, 1), c(1, 1, 1), list(c(-1, 1), c(1, 1, 1)))},
  
  and \code{list.snp} is set to
  
  \code{list(c(6, 7), c(3, 9, 10), list(c(2, 5), c(1, 4, 8)))}
  
  such that the interaction
  
  \code{(SNP6 != 1) & (SNP7 == 1)}
  
  is assumed to be explanatory for the first level of the three-categorical response, the interaction
  
  \code{(SNP3 == 1) & (SNP9 == 1) & (SNP10 == 1)}
  
  is assumed to be explanatory for the second level, and the interactions
  
  \code{(SNP2 != 1) & (SNP5 == 1)}\ \ \ and
  
  \code{(SNP1 == 1) & (SNP4 == 1) & (SNP8 == 1)},
  
  are assumed to be explanatory for the third level.
  
  For each of the levels, the probability of having this level given that an observation exhibits one, two, ...
  of the interactions intended to be explanatory for that level is determined using the corresponding
  logistic regression model. Afterwards, the value of the response for an observation showing one, two, ...
  of the interactions explanatory for a particular level is randomly drawn using the above probability \eqn{p}
  for the particular level and \eqn{(1-p)/(n_{cat}-1)}{(1-p)/(n.cat - 1)} as probabilities for the other
  \eqn{(n_{cat}-1)}{(n.cat - 1)} levels. If an observation exhibits none of the explanatory interactions,
  its response value is randomly drawn using the probabilities \eqn{\exp{beta0}/(1+\exp{beta0})}.
  
}
\value{
  An object of class \code{simSNPcatResponse} consisting of
  \item{x}{a matrix with \code{n.obs} rows and \code{n.snp} columns containing the simulated SNP values.}
  \item{y}{a vector of length \code{n.obs} composed of the values of the response.}
  \item{models}{a character vector naming the level-wise logistic regression models.}
  \item{maf}{a vector of length \code{n.snp} composed of the minor allele frequencies.}
  \item{tab.explain}{a data frame summarizing the results of the simulation.}
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{simulateSNPs}}, \code{\link{simulateSNPglm}}}
\examples{\dontrun{
# The simulated data set described in Details.

sim1 <- simulateSNPcatResponse()
sim1

# Specifying the values of the response by the levels with
# the largest probability.

sim2 <- simulateSNPcatResponse(sample.y = FALSE)
sim2

# If ((SNP4 != 2) & (SNP3 == 1)), (SNP5 ==3), and
# ((SNP12 !=1) & (SNP9 == 3)) should be the three interactions
# (or variables) that are explanatory for the three levels
# of the response, list.ia and list.snp are specified as follows.

list.ia <- list(c(-2, 1), 3, c(-1,3))
list.snp <- list(c(4, 3), 5, c(12,9))

# The categorical response and a data set consisting of 
# 800 observations and 25 SNPs, where the minor allele
# frequency of each SNP is randomly drawn from a
# uniform distribution with minimum 0.1 and maximum 0.4,
# is then generated by

sim3 <- simulateSNPcatResponse(n.obs = 800, n.snp = 25,
  list.ia = list.ia, list.snp = list.snp, maf = c(0.1, 0.4))
sim3

}}

\keyword{datagen}
