% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lineqGPlikelihoods.R
\name{logLikAdditiveFun}
\alias{logLikAdditiveFun}
\title{Log-Likelihood of a Additive Gaussian Process.}
\usage{
logLikAdditiveFun(
  par = unlist(purrr::map(model$kernParam, "par")),
  model,
  parfixed = NULL,
  mcmc.opts = NULL,
  estim.varnoise = FALSE
)
}
\arguments{
\item{par}{the values of the covariance parameters.}

\item{model}{an object with \code{"lineqAGP"} S3 class.}

\item{parfixed}{not used.}

\item{mcmc.opts}{not used.}

\item{estim.varnoise}{If \code{true}, a noise variance is estimated.}
}
\value{
The value of the negative log-likelihood.
}
\description{
Compute the negative log-likelihood of an Additive Gaussian Process.
}
\references{
Rasmussen, C. E. and Williams, C. K. I. (2005),
"Gaussian Processes for Machine Learning (Adaptive Computation and Machine Learning)".
\emph{The MIT Press}.
\href{http://www.gaussianprocess.org/gpml/}{[link]}
}
\seealso{
\code{\link{logLikAdditiveGrad}}
}
\author{
A. F. Lopez-Lopera.
}
