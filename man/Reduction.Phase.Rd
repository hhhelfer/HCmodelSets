\name{Reduction.Phase}
\alias{Reduction.Phase}

\title{Reduction by successive traversal of hypercubes proposed by Cox, D. R. & Battey, H. S. (2017)}
\description{
This function traverses successively lower dimensional hypercubes, discarding variables according to the appropriate decision rules. It provides the number and indices of variables selected at each stage.
}
\usage{
Reduction.Phase(X,Y,family=gaussian,
                dmHC=NULL,vector.signif=NULL,seed.HC = NULL, Cox.Hazard = FALSE)
}

\arguments{
  \item{X}{Design matrix.}
  \item{Y}{Response vector.}
  \item{family}{A description of the error distribution and link function to be used in the model. For glm this can be a character string naming a family function, a family function or the result of a call to a family function. See \code{\link{family}} for more details.}
  \item{dmHC}{Dimension of the hypercube to be used in the first-stage reduction. This version supports dimensions 2,3,4 and 5. If not specified a sensible value is calculated and used.}
  \item{vector.signif}{Vector of decision rules to be used at each stage of the reduction. The first value makes reference to the decision rule for the highest dimensional hypercube and so on. If values are less than 1, this specifies a significance level of a test. All variables significant at this level in at least half the analyses in which they appear will be retained. If the value is 1 or 2, variables are retained if they are among the 1 or 2 most significant in at least half the analyses in which they appear. If unspecified a default rule is used.}
  \item{seed.HC}{Seed for randomization of the variable indices in the hypercube. If not provided, the variables are arranged according to their original order.}
  \item{Cox.Hazard}{If TRUE fits proportional hazards regression model. The family argument will be ignored if Cox.Hazard=TRUE.}
}

\value{
\item{Matrix.Selection}{The number of variables selected at each reduction of the hypercube.}
\item{List.Selection}{The indexes of the variables selected at each reduction of the hypercube.}
}

\references{
Cox, D. R. and Battey, H. S. (2017). Large numbers of explanatory variables, a semi-descriptive analysis. \emph{Proceedings of the National Academy of Sciences}, 114(32), 8592-8595.

Battey, H. S. & Cox, D. R. (2018). Large numbers of explanatory variables: a probabilistic assessment. \emph{Proceedings of the Royal Society of London, A.}, 474(2215), 20170631.
}
\author{
Hoeltgebaum, H. H.
}

\section{Acknowledgement}{The work was supported by the UK Engineering and Physical Sciences Research Council under grant number EP/P002757/1.}

\examples{
\donttest{
## Generates a random DGP
dgp = DGP(s=5, a=3, sigStrength=1, rho=0.9, n=100, intercept=5, noise=1,
          var=1, d=1000, DGP.seed = 2018)

#Reduction Phase using only the first 70 observations
outcome.Reduction.Phase =  Reduction.Phase(X=dgp$X[1:70,],Y=dgp$Y[1:70],
                                           family=gaussian, seed.HC = 1012)

# Not run, using vector.signif argument
# Fixing a decision rule of getting the 2 most significant in the first reduction
# and in the subsequent reduction, only those variables significant at 0.001 level
# outcome.Reduction.Phase =  Reduction.Phase(X=dgp$X[1:70,],Y=dgp$Y[1:70],
#                           vector.signif = c(2,0.001), family=gaussian, dmHC = 3)

}
\dontshow{
dgp = DGP(s=5, a=3, sigStrength=1, rho=0.9, n=20, intercept=5, noise=1,
          var=1, d=50, DGP.seed = 2019)

#Reduction Phase using only the first 70 observations
outcome.Reduction.Phase =  Reduction.Phase(X=dgp$X[1:10,],Y=dgp$Y[1:10],
                                           dmHC = 2, family=gaussian, seed.HC = 1093)
}

}
