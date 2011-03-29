\name{abf}
\alias{abf}

\title{
Approximate Bayes Factor
}
\description{
Computes the Approximate Bayes Factor proposed by Wakefield (2009) for test statistics \code{theta / sqrt(V)}
that under the null hypothesis are assumed to follow an asymptotic standard normal distribution.
}
\usage{
abf(theta, V, W, numerator = 0, pi1 = NA)
}

\arguments{
  \item{theta}{a vector of numeric values, e.g., the maximum likelihood estimates for the parameter of
    a logistic regression model computed by separately applying this simple logistic regression to several SNPs.
    It is thus assumed that under the null hypothesis \code{theta / sqrt(V)} is asymptotically standard normal distributed. 
}
  \item{V}{a vector of the same length as \code{theta} containing the variances of the estimates comprised by \code{theta}.
}
  \item{W}{the prior variance. Must be either a positive value or a vector of the same length as \code{theta} consisting of
    positive values.
}
  \item{numerator}{either 0 or 1, specifying whether the numerator of the approximate Bayes factor comprises the 
    probability for the null hypothesis or the probability for the alternative hypothesis. 
}
  \item{pi1}{either a numeric value between 0 and 1 specifying the prior probability of association or a vector of the 
    same length as \code{theta} specifying for each of the SNPs a prior probability that this SNP is associated with the response. 
    If specified, prior odds, posterior odds,
    and depending on \code{numerator} the Bayesian False Discovery Probability (\code{numerator = 0}) or the posterior
    probability of association (\code{numerator = 1}) are computed. If \code{NA}, only the approximate Bayes factors are
    returned.
}
}

\value{
If \code{pi1 = NA}, a vector of the same length as \code{theta} containing the values of the approximate Bayes factor.
If \code{pi1} is specified, a list consisting of
\item{ABF}{a numeric vector containing the values of the approximate Bayes factors,}
\item{priorOdds}{either a numeric value or a numeric vector comprising the prior odds of association (if \code{numerator = 1})
  or no association (if \code{numerator = 0}),}
\item{postOdds}{a numeric vector containing the posterior odds of association (if \code{numerator = 1}) or no association
  (if \code{numerator = 0}),}
and either
\item{BFDP}{a numeric vector containing the Bayesian False Discovery Probabilities for the SNPs (if \code{numerator = 0}),}
or
\item{PPA}{a numeric vector comprising the posterior probabilities of association (if \code{numerator = 1})l} 
}

\references{
Wakefield, J. (2007). A Bayesian Measure of Probability of False Discovery in Genetic Epidemiology Studies. \emph{American Journal
  of Human Genetics}, 81, 208-227.
}
\author{
Holger Schwender, \email{holger.schw@gmx.de}
}

\keyword{htest}
