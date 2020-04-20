\name{LymphomaData}
\alias{patient.data}
\docType{data}
\title{
Lymphoma patients data set.
}
\description{
Data set of lymphoma patients used in the study of Alizadeh et al. (2000) and also Simon et al. (2011).
}
\usage{data(LymphomaData)}
\format{
  \describe{
    \item{\code{patient.data}}{A list with survival times, staus and covariates from patients.}
  }
}
\value{
\item{x}{Covariates from patients.}
\item{time}{Survival times.}
\item{status}{Patient status.}
}
\references{
Alizadeh, A. A., et al. (2000). Distinct types of diffuse large B-cell lymphoma identified by gene expression profiling. Nature, 403(6769), p.503.

Simon, N., Friedman, J., Hastie, T., & Tibshirani, R. (2011). Regularization paths for Cox’s proportional hazards model via coordinate descent. Journal of statistical software, 39(5), 1.
}
\examples{
data(LymphomaData)
x <- t(patient.data$x)
y <- patient.data$time
}
\keyword{datasets}
