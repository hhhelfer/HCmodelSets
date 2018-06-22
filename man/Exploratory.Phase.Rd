\name{Exploratory.Phase}
\alias{Exploratory.Phase}

\title{Perform the Exploratory phase on the hypercube dimension reduction proposed by Cox, D. R. & Battey, H. S. (2017)}
\description{
This function performs the exploratory phase on the variables retained through the reduction phase, returning any significant squared and interaction terms.
}
\usage{
Exploratory.Phase(X, Y, list.reduction, family=gaussian, log.transf=FALSE,
                  signif=0.01, silent=TRUE, Cox.Hazard = FALSE)
}

\arguments{
  \item{X}{Design matrix.}
  \item{Y}{Response vector.}
  \item{list.reduction}{Indices of retained variables from the reduction phase.}
  \item{family}{A description of the error distribution and link function to be used in the model. For glm this can be a character string naming a family function, a family function or the result of a call to a family function. See \code{\link{family}} for more details.}
  \item{log.transf}{Apply a log transformation to the elements of the design matrix. Default is FALSE.}
  \item{signif}{Significance level for the assessment of squared and interaction terms. The default is 0.01.}
  \item{silent}{By default, silent=TRUE. If silent=FALSE the user can decide upon the exclusion of individual interaction terms.}
  \item{Cox.Hazard}{If TRUE fits proportional hazards regression model. The family argument will be ignored if Cox.Hazard=TRUE.}
}

\value{
\item{mat.select.SQ}{Indices of variables with significant squared terms.}
\item{mat.select.INTER}{Indices of the pairs of variables with significant interaction terms.}
}

\references{
Cox, D. R., & Battey, H. S. (2017). Large numbers of explanatory variables, a semi-descriptive analysis. \emph{Proceedings of the National Academy of Sciences}, 114(32), 8592-8595.

Battey, H. S. & Cox, D. R. (2018). Large numbers of explanatory variables: a probabilistic assessment. \emph{Proceedings of the Royal Society of London, A.}, to appear
}
\author{
Hoeltgebaum, H. H.
}

\section{Acknowledgement}{The work was supported by the UK Engineering and Physical Sciences Research Council under grant number EP/P002757/1}

\seealso{
\code{\link{Reduction.Phase}}
}
\examples{
## Generates a random DGP using the DGP from the HCmodelSelection package
dgp = DGP(s=5, a=3, sigStrength=1, rho=0.9, n=100, intercept=5, noise=1,
          var=1, d=1000, DGP.seed = 2019)

## Design matrix
dgp$X

## Response variable
dgp$Y

#Reduction Phase using only the first 70 observations
outcome.Reduction.Phase =  Reduction.Phase(X=dgp$X[1:70,],Y=dgp$Y[1:70],
                                           dmHC = 5, family=gaussian, seed.HC = 1093)

# Exploratory Phase using only the first 70 observations, choosing the variables which
# were selected at least two times in the third dimension reduction

idxs = outcome.Reduction.Phase$List.Selection$`Hypercube with dim 2`$numSelected1
outcome.Exploratory.Phase =  Exploratory.Phase(X=dgp$X[1:70,],Y=dgp$Y[1:70],
                                               list.reduction = idxs,
                                               family=gaussian, signif=0.01)



}

