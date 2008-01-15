\name{analyse.models}
\alias{analyse.models}
\title{ Summarize MCMC sample of Bayesian logic regression models}

\description{
 For an MCMC sample of Bayesian logic regression models obtained with \code{fblr}
 the distribution of the model size and the most common logic predictors with up
 to three binaries are reported.  
}
\usage{
analyse.models(file, size.freq = TRUE, moco = c(20, 10), int.freq = TRUE, 
                kmax = 10, int.level = 2, bin.names = NULL)
}

\arguments{
  \item{file}{character string naming file where MCMC output of \code{fblr} is 
    stored.}
  \item{size.freq}{determines whether distribution of model size is reported as
    frequencies (default) or proportions.}
  \item{moco}{a vector of length 2 or 3 that determines how many of the most 
  common main effects, two-factor interactions and (possibly) three-factor 
  interactions are reported.}
  \item{int.freq}{determines whether the number (default) or the proportion of 
  models containing a specific interaction is reported.}
  \item{kmax}{the maximum number of allowed logic predictors used in \code{fblr}.}
  \item{int.level}{the maximum number of allowed binaries in a logic predictor 
  used in \code{fblr}.}
  \item{bin.names}{character vector of names for the binary variables. If no 
  names are supplied, binaries are referred to with their indices.}
}

\details{
The logic regression models visited during the MCMC run are stored by 
\code{fblr} in the rows of a matrix in the following fashion: 
Position 1 contains the number of logic predictors in the model. The next 
\code{kmax * (int.level + 1)} positions contain the predictors, each predictor being 
coded as c(number of binaries in predictor, indices of binaries), where negative
indices denote the complement of a variable. It follow the log-likelihood of 
the model, the value of the precision of the regression parameters and the 
\code{kmax+1} regression parameters. Zeros indicate empty entries. 
\code{analyse.models} extracts some of the most interesting information, namely
which logic predictors occur most often in the visited models, from the sample. 
The complement of a binary is indicated with a minus sign preceding its name.       
}
\value{
  \item{size}{table of model sizes.}
  \item{ones}{table of the \code{moco[1]} most common single-binary predictors.}
  \item{twos}{table of the \code{moco[2]} most common two-binaries predictors.}
  \item{threes}{table of the \code{moco[3]} most common three-binaries predictors.}
}
\author{ Arno Fritsch, \email{arno.fritsch@uni-dortmund.de}}
\seealso{\code{\link{fblr}}, \code{\link{predictFBLR}} }
\examples{\dontrun{
# Use fblr on some simulated SNP data
snp <- matrix(rbinom(500*20,2,0.3),ncol=20)
bin <- snp2bin(snp)
int <- apply(bin,1,function(x) (x[1] == 1 & x[3] == 0)*1)
case.prob <- exp(-0.5+log(5)*int)/(1+exp(-0.5+log(5)*int))
y <- rbinom(nrow(snp),1,prob=case.prob)
fblr(y, bin, niter=1000, nburn=0)

analyse.models("fblr_mcmc.txt")

# with SNP names
name.snp <- LETTERS[1:20]
name.bin <- paste(rep(name.snp,each=2), c("_d","_r"),sep="")

analyse.models("fblr_mcmc.txt", bin.names=name.bin)
  }
}
\keyword{regression}
\keyword{nonlinear}
