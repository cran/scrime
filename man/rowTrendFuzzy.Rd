\name{rowTrendFuzzy}
\alias{rowTrendFuzzy}
\alias{getMatFuzzy}

\title{
Trend Test for Fuzzy Genotype Calls
}
\description{
\code{rowTrendFuzzy} performs the trend test proposed by Louis et al. (2010) based on fuzzy genotype calls, i.e. the weighted sums
over the confidences for the three genotypes as they are determined by preprocessing algorithms (e.g., CRLMM)
or imputation procedures.

Given the confidences and scores for all three genotypes, \code{getMatFuzzy} constructs a matrix containing the fuzzy genotype calls.
}
\usage{
rowTrendFuzzy(score, probs, y, mat.fuzzy = NULL, 
    alternative = c("two.sided", "less", "greater"), 
    check = TRUE)
    
getMatFuzzy(score, probs, check = TRUE)
}
\arguments{
  \item{score}{either a numeric vector of length 2 or 3, or a character string. 
  
    If the latter, \code{score} must
    be either \code{"additive"}, \code{"dominant"}, \code{"recessive"}, or an abbreviation of these terms. If 
    \code{score = "additive"}, then \code{score} is set to \code{c(0, 1, 2)}. If \code{score = "dominant"},
    then \code{score} is set to \code{c(0, 1, 1)}. And if \code{score = "recessive"}, \code{score} is set to \code{c(0, 0, 1)}.
    
    If \code{score} is a numeric vector of length 2, then the first value must be the score for the heterozygous genotype,
    and the second value the score for the homozygous variant genotype. If \code{score} is of length, the first entry of
    this vector must be a zero (the score for the homozygous reference genotype), followed by the scores for the heterozygous
    and homozygous variant genotype.
}

  \item{probs}{a list of length 2 or 3 consisting of matrices of the same size. Each matrix must contain the confidences
    for one of the three genotypes, where each row in the matrix represents a SNP and each column a sample (which must be in the same
    order in all matrices). The matrices in \code{probs} correspond to the scores in \code{score}. Thus, if \code{probs} has
    length 2, then the first matrix must contain the confidences for the heterozygous genotype, and the second matrix the confidences
    for the homozygous variant genotype. All elementwise sums of these two matrices must be smaller than or equal to 1. 
    If \code{probs} has length 3, the first object must be a matrix containing the confidences
    for the homozygous reference genotype, followed by the two other matrices comprising the confidences for the heterozygous and
    homozygous variant genotype. All elementwise sums of the three matrices must be equal to 1. 
}
  \item{y}{a vector of zeros and ones specifying which of the samples in the matrices in \code{probs} are cases (\code{1}) and which
    are controls (\code{0}).
}
  \item{mat.fuzzy}{a matrix containing the fuzzy genotype calls. If specified, \code{score} and \code{probs} are not allowed to be
    specified in \code{rowTrendFuzzy}. 
    If \code{NULL}, \code{mat.fuzzy} is determined by employing \code{getMatFuzzy}, i.e.\ by multiplying the confidences 
    in \code{probs} with the corresponding
    scores in \code{score} and computing the elementwise sums over the resulting matrices with the weighted confidences.
}
  \item{alternative}{a character string specifying the alternative hypothesis. Must be one of \code{"two.sided"} (default), 
    \code{"greater"}, or \code{"less"}. Abbreviations (e.g.\ the initial letter) for these choices are also allowed.
}
  \item{check}{logical specifying whether the specified objects should be extensively checked. If \code{FALSE}, only some basic
    checks are done. It is highly recommended to use \code{check = TRUE}, although the checking takes much more computing time
    than the determination of the values of the trend test statistic.
}
}


\value{
For \code{getMatFuzzy}, a matrix containing the fuzzy genotype calls. For \code{rowTrendFuzzy}, a list consisting of 
\item{stat}{a vector containing the values of the trend test statistic for all SNPs comprised by \code{probs},}
\item{rawp}{a vector containing the unadjusted p-values computed for the values in \code{stat} based on a standard normal distribution
  and the specification of \code{alternative}.}
\item{theta}{a vector containing estimates for the log odds ratios for risk corresponding to \code{stat},}
\item{varTheta}{a vector containing the variance estimates for \code{theta}.} 
}

\references{
Louis, T.A., Carvalho, B.S., Fallin, M.D., Irizarry, R.A., Li, Q., and Ruczinski, I. (2010). Association Tests that Accommodate 
Genotyping Errors. In Bernardo, J.M., Bayarri, M.J., Berger, J.O., Dawid, A.P., Heckerman, D., Smith, A.F.M., and West, M. (eds.),
\emph{Bayesian Statistics 9}, 393-420. Oxford University Press, Oxford, UK. With Discussion. 
}


\seealso{
\code{\link{rowTrendStats}}, \code{\link{rowCATTs}}
}

\keyword{htest}
\keyword{array}
\keyword{manip}
