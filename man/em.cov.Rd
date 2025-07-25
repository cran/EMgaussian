% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emcov.R
\name{em.cov}
\alias{em.cov}
\title{EM algorithm for multivariate normal, covariance matrix parameterization}
\usage{
em.cov(
  dat,
  max.iter = 500,
  tol = 1e-05,
  start = c("diag", "pairwise", "listwise", "full"),
  debug = 0,
  ...
)
}
\arguments{
\item{dat}{Data frame or matrix that contains the raw data.}

\item{max.iter}{Max number of EM cycles.}

\item{tol}{Tolerance for change in parameter estimates across EM Cycles. If
all changes are less than \code{tol}, the algorithm terminates.}

\item{start}{Starting value method (see details).}

\item{debug}{(Experimental) set an integer value > 1 for some information as the algorithm runs.}

\item{...}{Space for additional arguments, not currently used.}
}
\value{
A list with the following:
\itemize{
 \item{\code{p.est}: all parameter estimates as a vector (means followed by
   unique elements of covariance matrix; each row below diagonal stacked).}
 \item{\code{mu}: estimated means.}
 \item{\code{S}: estimated covariance matrix.}
 \item{\code{it}: number of EM cycles completed.}
 \item{\code{conv}: boolean value indicating convergence (TRUE) or not (FALSE).}
}
}
\description{
EM algorithm for multivariate normal, covariance matrix parameterization
}
\details{
This function computes all means and covariances among a set of variables using
the Expectation-Maximization (EM) algorithm to handle missing values, and assuming
multivariate normality. The EM code was originally developed for the precision
matrix parameterization (\code{\link{em.prec}}), i.e., the parameters are the
means and the inverse of the covariance matrix. But, this is easily modifiable
to handle a covariance matrix parameterization such that means and covariances
are the model parameters.

For starting values, the function accepts either a list that has \code{mustart} and
\code{covstart} slots corresponding to the starting mean and covariance matrix. This
is useful if the user would like to use custom starting values. Otherwise, a
character corresponding to any of the options available in the
\code{\link{startvals.cov}} function will be used to take a guess at starting values.
}
\examples{
\donttest{
  library(psych)
  data(bfi)
  test <- em.cov(bfi[,1:25])
}

}
