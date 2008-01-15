\name{predictFBLR}
\alias{predictFBLR}
\title{ Predict Case Probabilities with Full Bayesian Logic Regression }
\description{
  Predicts case probabilities for binary data (usually SNP data dichotomized with
  \code{snp2bin}) based on an MCMC sample of Bayesian logic regression models
  obtained with \code{fblr}.
  }
\usage{
predictFBLR(file, bin, kmax = 10, int.level = 2)
}
\arguments{
  \item{file}{character string naming file where MCMC sample is stored.}
  \item{bin}{matrix of binary variables to make predictions for. One row is one
    observation. The number of binary variables has to be the same as used in \code{fblr}.}
  \item{kmax}{the maximum number of allowed logic predictors used in \code{fblr}.}
  \item{int.level}{the maximum number of allowed binaries in a logic predictor 
    used in \code{fblr}.}}
\value{
  Vector of length \code{nrow(bin)} with predicted case probabilities.
}
\author{ Arno Fritsch, \email{arno.fritsch@uni-dortmund.de}}
\seealso{\code{\link{fblr}}}
\examples{\dontrun{
# Use fblr on some simulated SNP data
snp <- matrix(rbinom(500 * 20, 2, 0.3), ncol = 20)
bin <- snp2bin(snp)
int <- apply(bin,1,function(x) (x[1] == 1 & x[3] == 0)*1)
case.prob <- exp(-0.5+log(5)*int)/(1+exp(-0.5+log(5)*int))
y <- rbinom(nrow(snp),1,prob=case.prob)
fblr(y, bin, niter=1000, nburn=0)

# Prediction for some new observations 
newbin <- snp2bin(matrix(rbinom(100 * 20, 2, 0.3), ncol = 20))
predictFBLR("fblr_mcmc.txt",newbin)

}
}
\keyword{regression}
\keyword{nonlinear}
