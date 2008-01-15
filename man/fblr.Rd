\name{fblr}
\alias{fblr}
\alias{fblr.weight}
\title{Full Bayesian Logic Regression for SNP Data}
\description{
  Performs full Bayesian logic regression for Single Nucleotide Polymorphism (SNP)
  data as described in Fritsch and Ickstadt (2007). 
  
  \code{fblr.weight} allows to incorporate prior pathway information by restricting 
  search for interactions to specific groups of SNPs and/or giving them different 
  weights. \code{fblr.weight} is only implemented for an interaction level of 2.
}
\usage{
fblr(y, bin, niter, thin = 5, nburn = 10000, int.level = 2, kmax = 10, 
  geo = 1, delta1 = 0.001, delta2 = 0.1, predict = FALSE, 
  file = "fblr_mcmc.txt")

fblr.weight(y, bin, niter, thin = 5, nburn = 10000, kmax = 10, geo = 1, 
  delta1 = 0.001, delta2 = 0.1, predict = FALSE, group = NULL, 
  weight = NULL, file = "fblr_mcmc.txt")
}
\arguments{
  \item{y}{binary vector indicating case-control status.}
  \item{bin}{binary matrix with number of rows equal to \code{length(y)}. 
  Usually the result of applying \code{snp2bin} to a matrix of SNP data.}
  \item{niter}{number of MCMC iterations after burn-in.}
  \item{thin}{after burn-in only every \code{thin}th iteration is kept.}
  \item{nburn}{number of burn-in iterations.}
  \item{int.level}{maximum number of binaries allowed in a logic predictor.
  Is fixed to 2 for \code{fblr.weight}.}
  \item{kmax}{maximum number of logic predictors allowed in the model.}
  \item{geo}{geometric penalty parameter for the number of binaries in a predictor.
  Value between 0 and 1. Default is \code{1}, meaning no penalty.}
  \item{delta1}{shape parameter for hierarchical gamma prior on precision of 
  regression parameters.}
  \item{delta2}{rate parameter for hierarchical gamma prior on precision of 
  regression parameters.}
  \item{predict}{should predicted case probabilities be returned?}
  \item{file}{character string naming a file to write the MCMC output to. If 
    \code{fblr} is called again, the file is overwritten.}
  \item{group}{list containing vectors of indices of binaries that are allowed 
  to interact. Groups may be overlapping, but every binary has to be in at least
  one group. Groups have to contain at least two binaries. Defaults to NULL,
  meaning that all interactions are allowed.}
  \item{weight}{vector of length \code{ncol(bin)} containing different relative 
  prior weights for binaries. Defaults to \code{NULL}, meaning equal weight for all binaries.}
}
\details{
The MCMC output in \code{file} can be analysed using the function 
\code{analyse.models}. In the help of this function it is also described how
the models are stored in \code{file}. 
}

\value{
  \item{accept}{acceptance rate of MCMC algorithm.}
  \item{pred}{vector of predicted case probabilities. Only given if 
  \code{predict = TRUE}.}

}
\references{ Fritsch, A. and Ickstadt, K.\ (2007). Comparing logic regression 
based methods for identifying SNP interactions. 
In \emph{Bioinformatics in Research and Development}, Hochreiter, S.\ and 
Wagner, R.\ (Eds.), Springer, Berlin.}

\author{Arno Fritsch, \email{arno.fritsch@uni-dortmund.de}}

\seealso{ \code{\link{analyse.models}},\code{\link{predictFBLR}}}
\examples{\dontrun{
# SNP dataset with 500 persons and 20 SNPs each,
# a two-SNP interaction influences the case probability
snp <- matrix(rbinom(500*20,2,0.3),ncol=20)
bin <- snp2bin(snp)
int <- apply(bin,1,function(x) (x[1] == 1 & x[3] == 0)*1)
case.prob <- exp(-0.5+log(5)*int)/(1+exp(-0.5+log(5)*int))
y <- rbinom(nrow(snp),1,prob=case.prob)

# normally more iterations should be used
fblr(y, bin, niter=1000, nburn=0)
analyse.models("fblr_mcmc.txt")

# Prior information: SNPs 1-10 belong to genes in one pathway, 
# SNPs 8-20 to another. Only interactions within a pathway are 
# considered and the first pathway is deemed to be twice as 
# important than the second.
fblr.weight(y,bin,niter=1000, nburn=0, group=list(1:20, 15:40), 
  weight=c(rep(2,20),rep(1,20)))
analyse.models("fblr_mcmc.txt")

}
}
\keyword{regression}
\keyword{nonlinear}
