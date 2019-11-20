#' @title  \code{"tmvrnorm"} Sampler for \code{"RSM"} (Rejection Sampling from the Mode) S3 Class
#' @description Sampler for truncated multivariate normal distributions
#' via RSM according to (Maatouk and Bay, 2017).
#' @param object an object with \code{"RSM"} S3 class containing:
#'        \code{mu} (mean vector), \code{Sigma} (covariance matrix),
#'        \code{lb} (lower bound vector), \code{ub} (upper bound vector).
#' @param nsim an integer corresponding to the number of simulations.
#' @param control extra parameters required for the MC/MCMC sampler.
#' @param ... further arguments passed to or from other methods.
#' @return A matrix with the simulated samples. Samples are indexed by columns.
#'
#' @seealso \code{\link{tmvrnorm.Gibbs}},
#'          \code{\link{tmvrnorm.HMC}}, \code{\link{tmvrnorm.ExpT}}
#' @author A. F. Lopez-Lopera.
#'
#' @references Maatouk, H. and Bay, X. (2017),
#' "Gaussian process emulators for computer experiments with inequality constraints".
#' \emph{Mathematical Geosciences},
#' 49(5):557-582.
#' \href{https://link.springer.com/article/10.1007/s11004-017-9673-2}{[link]}
#'
#' @examples
#' n <- 100
#' x <- seq(0, 1, length = n)
#' Sigma <- kernCompute(x1 = x, type = "gaussian", par = c(1,0.2))
#' tmgPar <- list(mu = rep(0,n), Sigma = Sigma + 1e-9*diag(n), lb = rep(-1,n), ub = rep(1,n))
#' class(tmgPar) <- "RSM"
#' y <- tmvrnorm(tmgPar, nsim = 10)
#' matplot(x, y, type = 'l', ylim = c(-1,1),
#'         main = "Constrained samples using RSM")
#' abline(h = c(-1,1), lty = 2)
#'
#' @importFrom MASS mvrnorm
#' @export
tmvrnorm.RSM <- function(object, nsim, control = NULL, ...) {
  # precomputing some terms according to Maatouk et al. [2016]
  tmvPar <- object
  if (!("map" %in% names(tmvPar)))
    tmvPar$map <- tmvPar$mu
  
  invSigma <- chol2inv(chol(tmvPar$Sigma))
  invSigmamu_star <-  invSigma %*% tmvPar$map

  # generating the samples
  xi <- matrix(0, nrow = length(tmvPar$map), ncol = nsim)
  for (i in seq(nsim)) {
    condition <- FALSE
    while (condition == FALSE) {
      xi_temp <- mvrnorm(n = 1, tmvPar$map, Sigma = tmvPar$Sigma)
      while(!all(xi_temp >= tmvPar$lb & xi_temp <= tmvPar$ub)) {
        xi_temp <- mvrnorm(n = 1, tmvPar$map, Sigma = tmvPar$Sigma)
      }
      t <- exp( t(tmvPar$map - xi_temp) %*% invSigmamu_star)
      u <- runif(1)
      if (u <= t) condition <- TRUE
    }
    xi[, i] <- xi_temp
  }
  return(xi)
}

# #' @title  \code{"tmvrnorm"} Sampler for \code{"RSMsov"} (Rejection Sampling from the Mode via SOV) S3 Class
# #' @description Sampler for truncated multivariate normal distributions
# #' via RSM according to (Maatouk and Bay, 2017).
# #' @param object an object with \code{"RSMsov"} S3 class containing:
# #'        \code{mu} (mean vector), \code{Sigma} (covariance matrix),
# #'        \code{lb} (lower bound vector), \code{ub} (upper bound vector).
# #' @param nsim an integer corresponding to the number of simulations.
# #' @param control extra parameters required for the MC/MCMC sampler.
# #' @param ... further arguments passed to or from other methods.
# #' @return A matrix with the simulated samples. Samples are indexed by columns.
# #'
# #' @seealso \code{\link{tmvrnorm.Gibbs}},
# #'          \code{\link{tmvrnorm.HMC}}, \code{\link{tmvrnorm.ExpT}}
# #' @author A. F. Lopez-Lopera.
# #'
# #' @references Maatouk, H. and Bay, X. (2017),
# #' "Gaussian process emulators for computer experiments with inequality constraints".
# #' \emph{Mathematical Geosciences},
# #' 49(5):557-582.
# #' \href{https://link.springer.com/article/10.1007/s11004-017-9673-2}{[link]}
# #'
# #' @examples
# #' n <- 100
# #' x <- seq(0, 1, length = n)
# #' Sigma <- kernCompute(x1 = x, type = "gaussian", par = c(1,0.2))
# #' tmgPar <- list(mu = rep(0,n), Sigma = Sigma + 1e-9*diag(n), lb = rep(-1,n), ub = rep(1,n))
# #' class(tmgPar) <- "RSMsov"
# #' y <- tmvrnorm(tmgPar, nsim = 10)
# #' matplot(x, y, type = 'l', ylim = c(-1,1),
# #'         main = "Constrained samples using RSM")
# #' abline(h = c(-1,1), lty = 2)
# #'
# #' @importFrom MASS mvrnorm
# #' @export
# tmvrnorm.RSMsov <- function(object, nsim, control = NULL, ...) {
#   # precomputing some terms according to Maatouk et al. [2016]
#   tmvPar <- object
#   if (!("map" %in% names(tmvPar)))
#     tmvPar$map <- tmvPar$mu
#   
#   Sigma <- tmvPar$Sigma
#   invSigma <- chol2inv(chol(tmvPar$Sigma))
#   invSigmamu_star <-  invSigma %*% tmvPar$map
#   
#   # condSd <- condVar <- rep(0, nrow(tmvPar$Sigma))
#   # condSd[1] <- sqrt(Sigma[1,1])
#   # for (k in 2:length(condSd)) {
#   #   Sigma_k <- Sigma[1:(k-1), 1:(k-1)]
#   #   invSigma_k <- chol2inv(chol(Sigma_k + 1e-9*diag(k-1)))
#   #   condVar[k] <- Sigma[k,k] -
#   #     Sigma[k,1:(k-1)] %*% invSigma_k %*% Sigma[1:(k-1),k]
#   #   condSd[k] <- sqrt(condVar[k])
#   # }
#   
#   mo <- tmvPar$map
#   A <- diag(length(tmvPar$map))
#   lb <- tmvPar$lb - A %*% mo
#   ub <- tmvPar$ub - A %*% mo
#   
#   Gamma <- A %*% Sigma %*% t(A)
#   Gamma <- Gamma + 1e-9*diag(nrow(Gamma))
#   d <- length(lb)
#   out <- cholperm(Gamma, lb, ub)
#   Lfull <- out$L
#   D <- diag(Lfull)
#   u <- out$u/D
#   l <- out$l/D
#   L <- out$L/D - diag(d)
#   perm <- out$perm
#   xmu <- nleq(l, u, L)
#   x <- xmu[1:(d - 1)]
#   nuexpt <- xmu[d:(2*d - 2)]
#   # psistar <- psy(x, L, l, u, nuexpt)
#   
#   # generating the samples
#   xi <- matrix(0, nrow = length(tmvPar$map), ncol = nsim)
#   for (i in seq(nsim)) {
#     condition <- FALSE
#     while (condition == FALSE) {
#       # xi_temp <- rep(0, length(tmvPar$map))
#       # for (k in 1:length(xi_temp)) {
#       #   u <- runif(1)
#       #   xi_temp[k] <-  tmvPar$map[k] + condSd[k]*qnorm(pnorm(tmvPar$lb[k], tmvPar$map[k], condSd[k]) +
#       #                         u*(pnorm(tmvPar$ub[k], tmvPar$map[k], condSd[k]) -
#       #                              pnorm(tmvPar$lb[k], tmvPar$map[k], condSd[k])))
#       # }
#       out <- mvnrnd(1, L, l, u, nuexpt)
#       logpr <- out$logpr
#       yexpt <- out$Z
#       # idx <- -log(runif(1)) > (psistar - logpr)
#       # sampexpt <-  out$Z[ , idx]
#       out <- sort(perm, decreasing = FALSE, index.return = TRUE)
#       order <- out$ix
#       yexpt <- yexpt[, 1:1]
#       yexpt <- Lfull %*% yexpt
#       yexpt <- yexpt[order, ]
#       # muexpt <- c(mo + qr.solve(A, c(nuexpt, nuexpt)))
#       xi_temp <- c(matrix(mo, nrow = nrow(Sigma), ncol = 1) + qr.solve(A, yexpt))
#       # # xi_temp <- mvrnorm(n = 1, tmvPar$map, Sigma = tmvPar$Sigma)
#       # # while(!all(xi_temp >= tmvPar$lb & xi_temp <= tmvPar$ub)) {
#       # #   xi_temp <- mvrnorm(n = 1, tmvPar$map, Sigma = tmvPar$Sigma)
#       # # }
#       # t <- exp( t(tmvPar$map - xi_temp) %*% invSigmamu_star)
#       # u2 <- runif(1)
#       # print(u2 <= t)
#       # if (u2 <= t) condition <- TRUE
#       condition <- TRUE
#     }
#     xi[, i] <- xi_temp
#   }
#   return(xi)
# }

#' @title  \code{"tmvrnorm"} Sampler for \code{"Gibbs"} (Gibbs Sampling) S3 Class
#' @description Sampler for truncated multivariate normal distributions
#' via Gibbs sampling using the package \code{restrictedMVN} (Taylor and Benjamini, 2017).
#' @param object an object with \code{"Gibbs"} S3 class containing:
#'        \code{mu} (mean vector), \code{Sigma} (covariance matrix),
#'        \code{lb} (lower bound vector), \code{ub} (upper bound vector).
#' @param nsim an integer corresponding to the number of simulations.
#' @param control extra parameters required for the MC/MCMC sampler.
#' @param ... further arguments passed to or from other methods.
#' @return A matrix with the simulated samples. Samples are indexed by columns.
#'
#' @seealso \code{\link{tmvrnorm.RSM}},
#'          \code{\link{tmvrnorm.HMC}}, \code{\link{tmvrnorm.ExpT}}
#' @author A. F. Lopez-Lopera.
#'
#' @references Taylor, J. and Benjamini, Y. (2017),
#' "RestrictedMVN: multivariate normal restricted by affine constraints".
#'
#' @examples
#' n <- 100
#' x <- seq(0, 1, length = n)
#' Sigma <- kernCompute(x1 = x, type = "gaussian", par = c(1,0.2))
#' tmgPar <- list(mu = rep(0,n), Sigma = Sigma + 1e-9*diag(n), lb = rep(-1,n), ub = rep(1,n))
#' class(tmgPar) <- "Gibbs"
#' y <- tmvrnorm(tmgPar, nsim = 10)
#' matplot(x, y, type = 'l', ylim = c(-1,1),
#'         main = "Constrained samples using Gibbs sampling")
#' abline(h = c(-1,1), lty = 2)
#'
#' @import restrictedMVN
#' @export
tmvrnorm.Gibbs <- function(object, nsim,
                           control = list(thinning = 1e2, burn.in = 1e2), ...) {
  if (!("thinning" %in% names(control)))
    control$thinning <- 1e2
  if (!("burn.in" %in% names(control)))
    control$burn.in <- 1e2
  tmvPar <- object
  if (!("map" %in% names(tmvPar)))
    tmvPar$map <- tmvPar$mu

  # simulating samples using the package "restrictedMVN"
  constr <- thresh2constraints(nrow(tmvPar$Sigma), tmvPar$lb, tmvPar$ub)
  xi <- sample_from_constraints(constr$linear_part, constr$offset, tmvPar$mu,
                                tmvPar$Sigma, initial_point = tmvPar$map,
                                ndraw = control$thinning*nsim,
                                burnin = control$burn.in)
  idx <- seq(control$thinning, control$thinning*nsim, control$thinning)
  return(t(xi[idx, ]))
}

#' @title  \code{"tmvrnorm"} Sampler for \code{"HMC"} (Hamiltonian Monte Carlo) S3 Class
#' @description Sampler for truncated multivariate normal distributions
#' via Hamiltonian Monte Carlo using the package \code{tmg}  (Pakman and Paninski, 2014).
#' @param object an object with \code{"HMC"} S3 class containing:
#'        \code{mu} (mean vector), \code{Sigma} (covariance matrix),
#'        \code{lb} (lower bound vector), \code{ub} (upper bound vector).
#' @param nsim an integer corresponding to the number of simulations.
#' @param control extra parameters required for the MC/MCMC sampler.
#' @param ... further arguments passed to or from other methods.
#' @return A matrix with the simulated samples. Samples are indexed by columns.
#'
#' @seealso \code{\link{tmvrnorm.RSM}}, \code{\link{tmvrnorm.Gibbs}},
#'          \code{\link{tmvrnorm.ExpT}}
#' @author A. F. Lopez-Lopera.
#'
#' @references Pakman, A. and Paninski, L. (2014),
#' "Exact Hamiltonian Monte Carlo for truncated multivariate Gaussians".
#' \emph{Journal of Computational and Graphical Statistics},
#' 23(2):518-542.
#' \href{https://www.tandfonline.com/doi/abs/10.1080/10618600.2013.788448}{[link]}
#'
#' @examples
#' n <- 100
#' x <- seq(0, 1, length = n)
#' Sigma <- kernCompute(x1 = x, type = "gaussian", par = c(1,0.2))
#' tmgPar <- list(mu = rep(0,n), Sigma = Sigma + 1e-9*diag(n), lb = rep(-1,n), ub = rep(1,n))
#' class(tmgPar) <- "HMC"
#' y <- tmvrnorm(tmgPar, nsim = 10)
#' matplot(x, y, type = 'l', ylim = c(-1,1),
#'         main = "Constrained samples using Hamiltonian MC")
#' abline(h = c(-1,1), lty = 2)
#'
#' @import tmg
#' @importFrom Matrix bdiag
#' @export
tmvrnorm.HMC <- function(object, nsim, control = list(burn.in = 1e2), ...) {
  if (!("burn.in" %in% names(control)))
    control$burn.in <- 1e2
  if (!("mvec" %in% names(control)))
    control$mvec <- length(object$mu)
  if (!("constrType" %in% names(control)))
    control$constrType <- "boundedness"
  tmvPar <- object
  if (!("map" %in% names(tmvPar)))
    tmvPar$map <- tmvPar$mu

  # precomputing some terms
  M <- init_vec <- g <- numeric()
  idxInit <- 1
  idxEnd <- 0
  epsilon <- 1e-9
  for (i in seq(length(control$constrType))) {
    idxInit <- idxEnd + 1
    idxEnd <- idxEnd + control$mvec[i]
    lsys <- bounds2lineqSys(control$mvec[i], tmvPar$lb[idxInit:idxEnd],
                            tmvPar$ub[idxInit:idxEnd], A = diag(control$mvec[i]),
                            lineqSysType = 'oneside')
    M <- bdiag(M, lsys$M)
    g <- c(g, lsys$g)
    init_vec <- c(init_vec, pmin(pmax(tmvPar$map[idxInit:idxEnd],
                                      tmvPar$lb[idxInit:idxEnd]+epsilon),
                                 tmvPar$ub[idxInit:idxEnd]-epsilon))
  }
  M <- M[,-1]

  # simulating samples using the package "tmg"
  H <- solve(tmvPar$Sigma)
  xi <- rtmg(nsim, M = H, r = as.vector(t(tmvPar$mu) %*% H), initial = init_vec,
             f = as.matrix(M), g = as.vector(g), burn.in = control$burn.in)
  return(t(xi))
}

#' @title  \code{"tmvrnorm"} Sampler for \code{"ExpT"} (Exponential Tilting) S3 Class
#' @description Sampler for truncated multivariate normal distributions
#' via exponential tilting using the package \code{TruncatedNormal} (Botev, 2017).
#' @param object an object with \code{"ExpT"} S3 class containing:
#'        \code{mu} (mean vector), \code{Sigma} (covariance matrix),
#'        \code{lb} (lower bound vector), \code{ub} (upper bound vector).
#' @param nsim an integer corresponding to the number of simulations.
#' @param control extra parameters required for the MC/MCMC sampler.
#' @param ... further arguments passed to or from other methods.
#' @return A matrix with the simulated samples. Samples are indexed by columns.
#'
#' @seealso \code{\link{tmvrnorm.RSM}}, \code{\link{tmvrnorm.Gibbs}},
#'          \code{\link{tmvrnorm.HMC}}
#' @author A. F. Lopez-Lopera.
#'
#' @references Botev, Z. I. (2017),
#' "The normal law under linear restrictions: simulation and estimation via minimax tilting".
#' \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)},
#' 79(1):125-148.
#' \href{https://rss.onlinelibrary.wiley.com/doi/pdf/10.1111/rssb.12162}{[link]}
#'
#' @examples
#' n <- 100
#' x <- seq(0, 1, length = n)
#' Sigma <- kernCompute(x1 = x, type = "gaussian", par = c(1,0.2))
#' tmgPar <- list(mu = rep(0,n), Sigma = Sigma + 1e-9*diag(n), lb = rep(-1,n), ub = rep(1,n))
#' class(tmgPar) <- "ExpT"
#' y <- tmvrnorm(tmgPar, nsim = 10)
#' matplot(x, y, type = 'l', ylim = c(-1,1),
#'         main = "Constrained samples using expontial tilting")
#' abline(h = c(-1,1), lty = 2)
#'
#' @export
tmvrnorm.ExpT <- function(object, nsim, control = NULL, ...) {
  tmvPar <- object
  # simulating samples using the package "TruncatedNormal" (Truncated Multivariate Normal)
  xi <- matrix(tmvPar$mu, nrow = nrow(tmvPar$Sigma), ncol = nsim) +
    TruncatedNormal::mvrandn(tmvPar$lb-tmvPar$mu, tmvPar$ub-tmvPar$mu, tmvPar$Sigma, nsim)
  return(xi)
}
