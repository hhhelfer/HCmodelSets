\name{DGP}
\alias{DGP}

\title{Data generating process used by Battey, H. S. & Cox, D. R. (2018).}
\description{
This function generates realizations of random variables as described in the simple example of Battey, H. S. & Cox, D. R. (2018).
}
\usage{
DGP(s,a,sigStrength,rho,n,noise=NULL,var,d,intercept,type.response="N",DGP.seed=NULL,
    scale=NULL,shape=NULL,rate=NULL)

}

\arguments{
  \item{s}{Number of signal variables.}
  \item{a}{Number of noise variables correlated with signal variables.}
  \item{sigStrength}{Signal strength.}
  \item{rho}{Correlation among signal variables and noise variables correlated with signal variables.}
  \item{n}{Sample size.}
  \item{noise}{Variance of the observations around the true regression line.}
  \item{var}{Variance of the potential explanatory variables.}
  \item{d}{Number of potential explanatory variables.}
  \item{intercept}{Expected value of the response variable when all potential explanatory variables are at zero. It is only considered when type.response="N".}
  \item{type.response}{Generates gaussian ("N") or survival ("S") data from a proportional hazards model with Weibull baseline hazard.}
  \item{DGP.seed}{Seed for the random number generator.}
  \item{scale}{scale parameter of the proportional hazards model with Weibull baseline hazard.}
  \item{shape}{shape parameter of the proportional hazards model with Weibull baseline hazard.}
  \item{rate}{rate parameter of the exponential distribution of censoring times. If not provided, uncensored data are generated.}
}

\value{
\item{X}{The simulated design matrix.}
\item{Y}{The simulated response variable.}
\item{TRUE.idx}{Indexes of the variables in the true model.}
\item{status}{If type.response="S", provides the status from survival data.}
}
\references{
Battey, H. S. & Cox, D. R. (2018). Large numbers of explanatory variables: a probabilistic assessment. \emph{Proceedings of the Royal Society of London, A.}, 474(2215), 20170631.
}
\author{
Hoeltgebaum, H. H.
}

\section{Acknowledgement}{The work was supported by the UK Engineering and Physical Sciences Research Council under grant number EP/P002757/1.}

\examples{
## Generates DGP
\donttest{
## Generates a random DGP
dgp = DGP(s=5, a=3, sigStrength=1, rho=0.9, n=100, intercept=5, noise=1,
          var=1, d=1000, DGP.seed = 2018)
          }
\dontshow{
dgp = DGP(s=5, a=3, sigStrength=1, rho=0.9, n=20, intercept=5, noise=1,
          var=1, d=50, DGP.seed = 2019)

}
}
