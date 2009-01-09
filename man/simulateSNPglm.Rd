\name{simulateSNPglm}
\alias{simulateSNPglm}

\title{Simulation of SNP data}
\description{
  Simulates SNP data. Interactions of some of the simulated SNPs
  are then used to specify either a binary or a
  quantitative response by a logistic or linear regression model, respectively.
}
\usage{
simulateSNPglm(n.obs = 1000, n.snp = 50, list.ia = NULL, list.snp = NULL, 
   beta0 = -0.5, beta = 1.5, maf = 0.25, sample.y = TRUE, p.cutoff = 0.5, 
   err.fun = NULL, rand = NA, ...)
}

\arguments{
  \item{n.obs}{number of observations that should be generated.}
  \item{n.snp}{number of SNPs that should be generated.}
  \item{list.ia}{a list consisting of numeric vectors (or values) specifying the genotypes
     of the SNPs that should be explanatory for the response. Each of these vectors
     must be composed of some of the numbers -3, -2, -1, 1, 2, 3, where 1 denotes
     the homozygous reference genotype, 2 the heterozygous genotype, and 3 the
     homozygous variant genotype, and a minus before these numbers means that
     the corresponding SNP should be not of this genotype. If, e.g., one of the
     vectors is given by \code{c(1, -1, -3)} and the corresponding vector in
     \code{list.snp} is \code{c(5, 7, 8)}, then the corresponding interaction
     used in the regression model to specify the response is
      
     \code{(SNP5 == 1) & (SNP7 != 1) & (SNP8 != 3)}. 
     
     For more details, see Details. Must be specified if \code{list.snp} is specified. 
     If both \code{list.ia} and \code{list.snp} are \code{NULL}, then the interactions
     shown in the Details section are used.}
  \item{list.snp}{a list consisting of numeric vectors specifying the SNPs that compose
     the interactions used in the regression model. Each of these vectors must have the
     same length as the corresponding vector in \code{list.ia}, and  must consist
     of integers between 1 and \code{n.snp}, where the integer \eqn{i} corresponds to
     the \eqn{i}th column of the simulated SNP matrix. If \code{list.ia} is specified but
     not \code{list.snp}, then the first \eqn{n} SNPs are used to generate the interactions,
     where \eqn{n} is the total number of values in \code{list.ia}. For the case that both
     \code{list.ia} and \code{list.snp} are not specified, see Details.}
  \item{beta0}{a numeric value specifying the intercept of the regression model.}
  \item{beta}{a non-negative numeric value or vector of the same length as \code{list.ia}
    (i.e.\ one numeric value for each interaction) specifying the parameters in the
    regression model.}
  \item{maf}{either an integer, or a vector of length 2 or \code{n.snp} specifying the
    minor allele frequency. If an integer, all the SNPs will have the same
    minor allele frequency. If a vector of length \code{n.snp}, each SNP will have the minor
    allele frequency specified in the corresponding entry of \code{maf}. If length 2, then
    \code{maf} is interpreted as the range of the minor allele frequencies, and for each SNP,
    a minor allele frequency will be randomly drawn from a uniform distribution with
    the range given by \code{maf}.}
  \item{sample.y}{should the values of the response in the logistic regression model be randomly drawn
    using the probabilities of the respective observations for being a case? If \code{FALSE},
    then the response value of an observation is 1 if its probability for being a case is larger
    than \code{p.cutoff}, and otherwise the observation is classified as a control (i.e.\ the value
    of the response is 0). Ignored if \code{err.fun} is specified.}
  \item{p.cutoff}{a probability, i.e.\ a numeric value between 0 and 1, naming the cutoff for an
    observation to be called a case if \code{sample.y = FALSE}. For details, see \code{sample.y}.
    Ignored if \code{sample.y = TRUE} or \code{err.fun} is specified.}
  \item{err.fun}{a function for generating the error that is added to the linear model to determine
    the value of the (quantitative) response. If \code{NULL}, a logistic regression model is fitted.
    If specified, a linear model is fitted. Therefore, this argument is used to differ between the
    two types of models. The specified function must have as first argument the number of values
    that should be generated and as output a vector consisting of these values. Further arguments
    can also be specified because of \code{\dots} in \code{simulateSNPglm}.
    If, e.g., \code{err.fun = rnorm}, then \code{rnorm(n.obs)} will be called to generate
    \code{n.obs} observations from a standard normal function.}
  \item{rand}{a numeric value for setting the random number generator in a reproducible state.}
  \item{\dots}{further arguments of the function specified by \code{err.fun}.}
}
\details{
  \code{simulateSNPglm} first simulates a matrix consisting of \code{n.obs} observations and \code{n.snp}
  SNPs, where the minor allele frequencies of these SNPs are given by \code{maf}.
  
  Note that all SNPs are currently simulated independently of each other such that they are unlinked.
  
  Afterwards, the response is determined by a regression model using the specifications of \code{list.ia},
  \code{list.snp}, \code{beta0} and \code{beta}. Depending on whether \code{err.fun} is specified or not,
  a linear or a logistic regression model is used, respectively, i.e.\ the response \eqn{Y} is continuous or binary.
  
  By default, a logistic regression model
  
  logit(Prob(\eqn{Y = 1})) = \code{beta0} + \code{beta[1]} * \eqn{L_1}{L1} + \code{beta[2]} * \eqn{L_2}{L2} + \dots
  
  is fitted, since \code{err.fun = NULL}. 
  
  If both \code{list.ia} and \code{list.snp} are \code{NULL}, then interactions similar to the one
  considered, e.g., in Nunkesser et al.\ (2007) or Schwender et al.\ (2007) are used, i.e.\
  
  \eqn{L_1}{L1} = \code{(SNP6 != 1) & (SNP7 == 1)}\ \ \ and 
  
  \eqn{L_2}{L2} = \code{(SNP3 == 1) & (SNP9 == 1) & (SNP10 == 1)},
  
  by setting \code{list.ia = list(c(-1, 1), c(1, 1, 1))} and
  \code{list.snp = list(c(6, 7), c(3, 9, 10))}. 
  
  Using the above model Prob(\eqn{Y = 1}) is computed for each observation, and its value of the response is
  determined either by a random draw from a Bernoulli distribution using this probability (if \code{sample.y = TRUE}),
  or by evaluating if Prob(\eqn{Y = 1}) \eqn{>} \code{p.cutoff} (if \code{sample.y = FALSE}).
  
  If \code{err.fun} is specified, then the linear model
  
  \eqn{Y} = \code{beta0} + \code{beta[1]} * \eqn{L_1}{L1} + \code{beta[2]} * \eqn{L_2}{L2} + \dots + error
  
  is used to determine the values of the response \eqn{Y}, where the values for error are given
  by the output of a call of \code{err.fun}.
}
\value{
  An object of class \code{simSNPglm} consisting of
  \item{x}{a matrix with \code{n.obs} rows and \code{n.snp} columns containing the simulated SNP values.}
  \item{y}{a vector of length \code{n.obs} composed of the values of the response.}
  \item{beta0}{the value of the intercept.}
  \item{beta}{the vector of parameters.}
  \item{ia}{a character vector naming the explanatory interactions.}
  \item{maf}{a vector of length \code{n.snp} composed of the minor allele frequencies.}
  \item{prob}{a vector of length \code{n.obs} consisting of the values of Prob(\eqn{Y = 1}) (if \code{err.fun = NULL}).}
  \item{err}{a vector of length \code{n.obs} composed of the values of the error (if \code{err.fun} is specified).}
  \item{p.cutoff}{the value of \code{p.cutoff} (if \code{err.fun = NULL} and \code{sample.y = FALSE}).}
  \item{err.call}{a character string naming the call of the error function (if \code{err.fun} is specified).}
}
\references{
Nunkesser, R., Bernholt, T., Schwender, H., Ickstadt, K. and Wegener, I.\ (2007).
Detecting High-Order Interactions of Single Nucleotide Polymorphisms Using Genetic Programming.
\emph{Bioinformatics}, 23, 3280-3288.

Schwender, H.\ (2007). Statistical Analysis of Genotype and Gene Expression Data. \emph{Dissertation},
Department of Statistics, University of Dortmund. 
}

\author{Holger Schwender, \email{holger.schwender@udo.edu}}

\seealso{\code{\link{simulateSNPs}}, \code{\link{summary.simSNPglm}}, \code{\link{simulateSNPcatResponse}}}

\examples{\dontrun{
# The simulated data set described in Details.

sim1 <- simulateSNPglm()
sim1

# A bit more information: Table of probabilities of being a case
# vs. numbers of cases and controls.

summary(sim1)

# Calling an observation a case if its probability of being
# a case is larger than 0.5 (the default for p.cutoff).

sim2 <- simulateSNPglm(sample.y = FALSE)
summary(sim2)

# If ((SNP4 != 2) & (SNP3 == 1)), (SNP5 ==3) and
# ((SNP12 !=1) & (SNP9 == 3)) should be the three interactions
# (or variables) that are explanatory for the response,
# list.ia and list.snp are specified as follows.

list.ia <- list(c(-2, 1), 3, c(-1,3))
list.snp <- list(c(4, 3), 5, c(12,9))

# The binary response and the data set consisting of 
# 600 observations and 25 SNPs, where the minor allele
# frequency of each SNP is randomly drawn from a
# uniform distribution with minimum 0.1 and maximum 0.4,
# is then generated by

sim3 <- simulateSNPglm(n.obs = 600, n.snp = 25,
  list.ia = list.ia, list.snp = list.snp, maf = c(0.1, 0.4))
sim3

summary(sim3)
  
# If the response should be quantitative, err.fun has
# to be specified. To use a normal distribution with mean 0
# (default in rnorm) and a standard deviation of 2 
# as the distribution of the error, call

simulateSNPglm(err.fun = rnorm, sd = 2)
 

}}

\keyword{datagen}
