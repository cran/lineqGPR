% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lineqGPsamplers.R
\name{tmvrnorm.HMC}
\alias{tmvrnorm.HMC}
\title{\code{"tmvrnorm"} Sampler for \code{"HMC"} (Hamiltonian Monte Carlo) S3 Class}
\usage{
\method{tmvrnorm}{HMC}(object, nsim, control = list(burn.in = 100), ...)
}
\arguments{
\item{object}{an object with \code{"HMC"} S3 class containing:
\code{mu} (mean vector), \code{Sigma} (covariance matrix),
\code{lb} (lower bound vector), \code{ub} (upper bound vector).}

\item{nsim}{an integer corresponding to the number of simulations.}

\item{control}{extra parameters required for the MC/MCMC sampler.}

\item{...}{further arguments passed to or from other methods.}
}
\value{
A matrix with the simulated samples. Samples are indexed by columns.
}
\description{
Sampler for truncated multivariate normal distributions
via Hamiltonian Monte Carlo using the package \code{tmg}  (Pakman and Paninski, 2014).
}
\examples{
n <- 100
x <- seq(0, 1, length = n)
Sigma <- kernCompute(x1 = x, type = "gaussian", par = c(1,0.2))
tmgPar <- list(mu = rep(0,n), Sigma = Sigma + 1e-9*diag(n), lb = rep(-1,n), ub = rep(1,n))
class(tmgPar) <- "HMC"
y <- tmvrnorm(tmgPar, nsim = 10)
matplot(x, y, type = 'l', ylim = c(-1,1),
        main = "Constrained samples using Hamiltonian MC")
abline(h = c(-1,1), lty = 2)

}
\references{
Pakman, A. and Paninski, L. (2014),
"Exact Hamiltonian Monte Carlo for truncated multivariate Gaussians".
\emph{Journal of Computational and Graphical Statistics},
23(2):518-542.
\href{https://www.tandfonline.com/doi/abs/10.1080/10618600.2013.788448}{[link]}
}
\seealso{
\code{\link{tmvrnorm.RSM}}, \code{\link{tmvrnorm.ExpT}}
}
\author{
A. F. Lopez-Lopera.
}
