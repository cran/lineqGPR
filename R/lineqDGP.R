#' @title Kernel Matrix for \code{"lineqDGP"} Models.
#' @description Compute the kernel matrix for \code{"lineqDGP"} models.
#' @param u a discretization vector of the input locations.
#' @param constrType a character string corresponding to the type of the inequality constraint.
#' Options: "boundedness", "monotonicity", "convexity".
#' @param kernType a character string corresponding to the type of the kernel.
#' Options: "gaussian", "matern32", "matern52", "exponential".
#' @param par the values of the kernel parameters (variance, lengthscale).
#' @param d a number corresponding to the dimension of the input space.
#' @return Kernel matrix \eqn{K(u,u)}
#'
#' @seealso \code{\link{kernCompute}}
#' @author A. F. Lopez-Lopera.
#'
#' @references Maatouk, H. and Bay, X. (2017),
#' "Gaussian process emulators for computer experiments with inequality constraints".
#' \emph{Mathematical Geosciences},
#' 49(5):557-582.
#' \href{https://link.springer.com/article/10.1007/s11004-017-9673-2}{[link]}
#'
#' @examples
#' x <- seq(0, 1, 0.01)
#' par <- c(1, 0.1)
#' Kb <- ineqConstrKernCompute(x, constrType = "boundedness", kernType = "gaussian", par)
#' image(Kb, main = "covariance matrix for boundedness constraints")
#' Km <- ineqConstrKernCompute(x, constrType = "monotonicity", kernType = "gaussian", par)
#' image(Km, main = "covariance matrix for monotonicity constraints")
#' Kc <- ineqConstrKernCompute(x, constrType = "convexity", kernType = "gaussian", par)
#' image(Kc, , main = "covariance matrix for convexity constraints")
#'
#' @export
ineqConstrKernCompute <- function(u,
                                  constrType = c("boundedness","monotonicity","convexity"),
                                  kernType, par, d = 1L) {
  # defining the kernel function
  kernName <- paste("k", d, kernType, sep = "")
  kernFun <- try(get(kernName))
  if (class(kernFun) == "try-error") {
    stop('kernel "', kernType, '" is not supported')
  } else {
    constrType <- match.arg(constrType)
    switch(constrType,
           boundedness = {
             kern <- kernFun(u, u, par, d)
           }, monotonicity = {
             kern11 <- kernFun(0, 0, par, d)
             kern21 <- attr(kernFun(u, 0, par, d), "derivative")$x1
             kern12 <- attr(kernFun(0, u, par, d), "derivative")$x2
             kern22 <- attr(kernFun(u, u, par, d), "derivative")$x1x2
             kern <- rbind(cbind(kern11, kern12),
                           cbind(kern21, kern22))
           }, convexity = {
             kern11 <- kernFun(0, 0, par, d)
             kern21 <- attr(kernFun(0, 0, par, d), "derivative")$x1
             kern31 <- attr(kernFun(u, 0, par, d), "derivative")$x1x1
             kern12 <- attr(kernFun(0, 0, par, d), "derivative")$x2
             kern22 <- attr(kernFun(0, 0, par, d), "derivative")$x1x2
             kern32 <- attr(kernFun(u, 0, par, d), "derivative")$x1x1x2
             kern13 <- attr(kernFun(0, u, par, d), "derivative")$x2x2
             kern23 <- attr(kernFun(0, u, par, d), "derivative")$x2x2x1
             kern33 <- attr(kernFun(u, u, par, d), "derivative")$x1x1x2x2
             kern <- rbind(cbind(kern11, kern12, kern13),
                           cbind(kern21, kern22, kern23),
                           cbind(kern31, kern32, kern33))
           })
    return(kern)
  }
}

#' @title Basis Functions for \code{"lineqDGP"} Models
#' @description Evaluate the basis functions for \code{"lineqDGP"} models.
#' @param x a vector with the input data.
#' @param m the number of basis functions used in the approximation.
#' @param d a number corresponding to the dimension of the input space.
#' @param constrType a character string corresponding to the type of the inequality constraint.
#' Options: "boundedness", "monotonicity", "convexity".
#' @return A matrix with the basis functions. The basis functions are indexed by rows.
#'
#' @author A. F. Lopez-Lopera.
#'
#' @references Maatouk, H. and Bay, X. (2017),
#' "Gaussian process emulators for computer experiments with inequality constraints".
#' \emph{Mathematical Geosciences},
#' 49(5):557-582.
#' \href{https://link.springer.com/article/10.1007/s11004-017-9673-2}{[link]}
#'
#' @examples
#' x <- seq(0, 1, 0.01)
#' Phib <- basisCompute.lineqDGP(x, m = 5, constrType = "boundedness")
#' matplot(Phib, type = "l", lty = 2,
#'         main = "Basis functions for boundedness constraints")
#' Phim <- basisCompute.lineqDGP(x, m = 5, constrType = "monotonicity")
#' matplot(Phim, type = "l", lty = 2,
#'         main = "Basis functions for monotonicity constraints")
#' Phic <- basisCompute.lineqDGP(x, m = 5, constrType = "convexity")
#' matplot(Phic, type = "l", lty = 2,
#'         main = "Basis functions for convexity constraints")
#'
#' @export
basisCompute.lineqDGP <- function(x, m, d = 1,
                                  constrType = c("boundedness","monotonicity","convexity")) {
  # precomputing some term
  if (!is.matrix(x) || ncol(x) != d)
    x <- matrix(x, ncol = 1)
  n <- nrow(x)
  delta <- 1/(m-1)
  u <- as.matrix(seq(0, 1, by = delta)) # discretization vector

  # computing the basis functions according to the type of inequality constraint
  constrType <- match.arg(constrType)
  switch(constrType,
         boundedness = {
           distAbs <- abs(outer(x[,1]/delta, u[,1]/delta, "-"))
           idx = distAbs <= 1
           Phi <- matrix(0, n, m)
           Phi[idx] <- 1 - distAbs[idx]
         }, monotonicity = {
           dist <- outer(x[,1], u[,1], "-")
           dist2 <- outer(x[,1], u[,1], "-")^2
           idx1 <- dist >= -delta & dist < 0
           idx2 <- dist >= 0 & dist < delta
           idx3 <- dist >= delta
           Phi <- matrix(0, n, m)
           Phi.fixed <- 0.5*delta + dist
           Phi[idx1] <- Phi.fixed[idx1] + 0.5*dist2[idx1]/delta
           Phi[idx2] <- Phi.fixed[idx2] - 0.5*dist2[idx2]/delta
           Phi[idx3] <- delta
           Phi[,1] <- Phi[,1] - 0.5*delta
         }, convexity = {
           dist <- outer(x[,1], u[,1], "-")
           dist2 <- outer(x[,1], u[,1], "-")^2
           dist3 <- outer(x[,1], u[,1], "-")^3
           idx1 <- dist >= -delta & dist < 0
           idx2 <- dist >= 0 & dist < delta
           idx3 <- dist >= delta
           Phi <- matrix(0, n, m)
           Phi.fixed <- delta^2/6 + 0.5*delta*dist + 0.5*dist2
           Phi[idx1] <- Phi.fixed[idx1] + (1/6)*dist3[idx1]/delta
           Phi[idx2] <- Phi.fixed[idx2] - (1/6)*dist3[idx2]/delta
           Phi[idx3] <- delta*dist[idx3]
           Phi[,1] <- 0.5*(Phi[,1] - delta^2/6)
         })
  return(Phi)
}

#' @title Linear Systems of Inequalities for \code{"lineqDGP"} Models
#' @description Build the linear system of inequalities for \code{"lineqDGP"} models.
#' @param d the number of linear inequality constraints.
#' @param constrType a character string corresponding to the type of the inequality constraint.
#' Options: "boundedness", "monotonicity", "convexity".
#' @param l the value (or vector) with the lower bound.
#' @param u the value (or vector) with the upper bound.
#' @param lineqSysType a character string corresponding to the type of the
#' linear system. Options: \code{twosides}, \code{oneside} (see \code{\link{bounds2lineqSys}} for more details).
#' @param rmInf If \code{TRUE}, inactive constraints are removed
#' (e.g. \eqn{-\infty \leq x \leq \infty}{-Inf \le x \le Inf}).
#' @return  A list with the linear system of inequalities: \code{list(A,l,u)} (\code{twosides}) or \code{list(M,g)} (\code{oneside} ).
#'
#' @seealso \code{\link{bounds2lineqSys}}
#' @author A. F. Lopez-Lopera.
#'
#' @references Maatouk, H. and Bay, X. (2017),
#' "Gaussian process emulators for computer experiments with inequality constraints".
#' \emph{Mathematical Geosciences},
#' 49(5):557-582.
#' \href{https://link.springer.com/article/10.1007/s11004-017-9673-2}{[link]}
#'
#' @examples
#' linSys1 <- lineqDGPSys(d = 5, constrType = "boundedness", l = 0, u = 1, lineqSysType = "twosides")
#' linSys1
#' linSys2 <- lineqDGPSys(d = 5, constrType = "boundedness", l = 0, u = 1, lineqSysType = "oneside")
#' linSys2
#'
#' @export
lineqDGPSys <- function(d,
                        constrType = c("boundedness","monotonicity","convexity"),
                        l = -Inf, u = Inf, lineqSysType = "twosides",
                        rmInf = TRUE) {
  constrType <- match.arg(constrType)
  switch(constrType,
         boundedness = {
           linSys <- bounds2lineqSys(d, l, u, lineqSysType = lineqSysType, rmInf = rmInf)
         }, monotonicity = {
           if (length(l) == 1)
             l <- c(-Inf, rep(l, d-1))
           linSys <- bounds2lineqSys(d, l, u, lineqSysType = lineqSysType, rmInf = rmInf)
         }, convexity = {
           if (length(l) == 1)
             l <- c(-rep(Inf,2), rep(l, d-2))
           linSys <- bounds2lineqSys(d, l, u, lineqSysType = lineqSysType, rmInf = rmInf)
         })
  return(linSys)
}

#' @title Creation Method for the \code{"lineqDGP"} S3 Class
#' @description Creation method for the \code{"lineqDGP"} S3 class.
#' @param x a vector or matrix with the input data. The dimensions should be indexed by columns.
#' @param y a vector with the output data.
#' @param constrType a character string (or list) corresponding to the type(s)
#' of inequality constraint(s).
#' Options: "boundedness", "monotonicity", "convexity", "linear".
#' @return A list with the following elements.
#' \item{x,y,constrType}{see \bold{Arguments}.}
#' \item{d}{a number corresponding to the input dimension.}
#' \item{localParam}{a list with specific parameters required for \code{"lineqGP"} models:
#' \code{m} (number of basis functions), \code{sampler}, and \code{samplingParams}.
#' See \code{\link{simulate.lineqDGP}}.}
#' \item{kernParam}{a list with the kernel parameters:
#' \code{par} (kernel parameters), \code{type}, \code{nugget}.
#' See \code{\link{kernCompute}}}
#' \item{bounds}{the limit values if \code{constrType = "boundedness"}.}
#'
#' @seealso \code{\link{augment.lineqDGP}}, \code{\link{predict.lineqDGP}},
#'          \code{\link{simulate.lineqDGP}}
#' @author A. F. Lopez-Lopera.
#'
#' @references Maatouk, H. and Bay, X. (2017),
#' "Gaussian process emulators for computer experiments with inequality constraints".
#' \emph{Mathematical Geosciences},
#' 49(5):557-582.
#' \href{https://link.springer.com/article/10.1007/s11004-017-9673-2}{[link]}
#'
#' @examples
#' # creating the model
#' sigfun <- function(x) return(1/(1+exp(-7*(x-0.5))))
#' x <- seq(0, 1, length = 5)
#' y <- sigfun(x)
#' model <- create(class = "lineqDGP", x, y, constrType = "monotonicity")
#' model
#'
#' @method create lineqDGP
#' @export
create.lineqDGP <- function(x, y,
                            constrType = c("boundedness","monotonicity","convexity")) {
  # changing the data as matrices
  if (!is.matrix(x))
    x <- as.matrix(x)
  if (!is.matrix(y) || ncol(y) != 1)
    y <- matrix(y, ncol = 1)
  d <- ncol(x) # dimension of the input space

  # creating some lists for the model
  kernParam <- list(par = c(sigma2 = 1^2, theta = 0.2),
                    type = "gaussian", nugget = 1e-7*sd(y))
  localParam <- list(m = 30, sampler = "HMC",
                     samplingParams = c(thinning = 1, burn.in = 1, scale = 0.1))

  # creating the full list for the model
  constrType <- match.arg(constrType)
  model <- list(x = x, y = y, constrType = constrType,
                d = d, localParam = localParam, kernParam = kernParam)
  if (any(constrType == "boundedness"))
    model$bounds <- c(lower = min(y) - 0.05*abs(max(y) - min(y)),
                      upper = max(y) + 0.05*abs(max(y) - max(y)))
  return(model)
}

#' @title Augmenting Method for the \code{"lineqDGP"} S3 Class
#' @description Augmenting method for the \code{"lineqDGP"} S3 class.
#' @param x an object with class \code{lineqDGP}.
#' @param ... further arguments passed to or from other methods.
#' @return An expanded \code{"lineqDGP"} object with the following additional elements.
#' \item{Phi}{a matrix corresponding to the hat basis functions.
#' The basis functions are indexed by rows.}
#' \item{Gamma}{the covariance matrix of the Gassian vector \eqn{\boldsymbol{\xi}}{\xi}.}
#' \item{(Lambda,lb,ub)}{the linear system of inequalities.}
#' \item{...}{further parameters passed to or from other methods.}
#'
#' @seealso \code{\link{create.lineqDGP}}, \code{\link{predict.lineqDGP}}, \code{\link{simulate.lineqDGP}}
#' @author A. F. Lopez-Lopera.
#'
#' @references Maatouk, H. and Bay, X. (2017),
#' "Gaussian process emulators for computer experiments with inequality constraints".
#' \emph{Mathematical Geosciences},
#' 49(5):557-582.
#' \href{https://link.springer.com/article/10.1007/s11004-017-9673-2}{[link]}
#'
#' @examples
#' # creating the model
#' sigfun <- function(x) return(1/(1+exp(-7*(x-0.5))))
#' x <- seq(0, 1, length = 5)
#' y <- sigfun(x)
#' model <- create(class = "lineqDGP", x, y, constrType = "monotonicity")
#'
#' # updating and expanding the model
#' model$localParam$m <- 30
#' model$kernParam$par <- c(1, 0.2)
#' model2 <- augment(model)
#' image(model2$Gamma, main = "covariance matrix")
#'
#' @importFrom broom augment
#' @export
augment.lineqDGP <- function(x, ...) {
  model <- x
  if (!("nugget" %in% names(model$kernParam)))
    model$kernParam$nugget <- 0
  if ("bounds" %in% names(model)) {
    bounds <- model$bounds
  } else {
    bounds <- c(0, Inf)
  }
  # passing some terms from the model
  x <- model$x
  m <- model$localParam$m
  u <- matrix(seq(0, 1, by = 1/(m-1)), ncol = 1) # discretization vector

  # computing the kernel matrix for the prior
  Phi <- basisCompute.lineqDGP(x, m, d = 1, model$constrType)
  Gamma <- ineqConstrKernCompute(u, model$constrType, model$kernParam$type,
                                 model$kernParam$par, d = 1)
  model$Gamma <- Gamma + model$kernParam$nugget*diag(nrow(Gamma))

  # precomputing some term according to Maatouk et al. [2016]
  mt <- nrow(model$Gamma) # number of basis functions (m) + number of initial conditions
  switch(model$constrType,
         boundedness = {
           model$Phi <- Phi
         }, monotonicity = {
           model$Phi <- cbind(matrix(1, nrow = nrow(x)), Phi)
         }, convexity = {
           model$Phi <- cbind(matrix(1, nrow = nrow(x)), x, Phi)
         }, {
           stop('constraint "', model$constrType, '" is not supported')
         }
  )
  lsys <- lineqDGPSys(mt, model$constrType, bounds[1], bounds[2], lineqSysType = "oneside")
  lsys2 <- lineqDGPSys(mt, model$constrType, bounds[1], bounds[2], rmInf = FALSE)

  # adding the parameters to the model structure
  model$lb <- lsys2$l
  model$ub <- lsys2$u
  model$lineqSys$M <- lsys$M  # for QP solve
  model$lineqSys$g <- -matrix(lsys$g)  # for QP solve
  model$localParam$mvec <- nrow(lsys2$A) # for HMC sampler
  return(model)
}

#' @title Prediction Method for the \code{"lineqDGP"} S3 Class
#' @description Prediction method for the \code{"lineqDGP"} S3 class.
#' @param object an object with class \code{"lineqDGP"}.
#' @param xtest a vector (or matrix) with the test input design.
#' @param ... further arguments passed to or from other methods.
#' @return An object with the predictions of \code{"lineqDGP"} models.
#' \item{lb}{The lower bound vector of the inequalities constraints.}
#' \item{ub}{The upper bound vector of the inequalities constraints.}
#' \item{Phi.test}{A matrix corresponding to the hat basis functions evaluated at \code{xtest}.
#' The basis functions are indexed by rows.}
#' \item{mu}{The unconstrained GP mean predictor.}
#' \item{xi.map}{The GP maximum a posteriori (MAP) predictor given the inequality constraints.}
#' \item{Sigma.xi}{The unconstrained GP prediction conditional covariance matrix.}
#'
#' @details The posterior sample-path of the finite-dimensional GP with inequality constraints is computed
#' according to (Maatouk and Bay, 2017).
#' @seealso \code{\link{create.lineqDGP}}, \code{\link{augment.lineqDGP}}, \code{\link{simulate.lineqDGP}}
#' @author A. F. Lopez-Lopera.
#'
#' @references Maatouk, H. and Bay, X. (2017),
#' "Gaussian process emulators for computer experiments with inequality constraints".
#' \emph{Mathematical Geosciences},
#' 49(5):557-582.
#' \href{https://link.springer.com/article/10.1007/s11004-017-9673-2}{[link]}
#'
#' @examples
#' # creating the model
#' sigfun <- function(x) return(1/(1+exp(-7*(x-0.5))))
#' x <- seq(0, 1, length = 5)
#' y <- sigfun(x)
#' model <- create(class = "lineqDGP", x, y, constrType = "monotonicity")
#'
#' # updating and expanding the model
#' model$localParam$m <- 30
#' model$kernParam$par <- c(1, 0.2)
#' model <- augment(model)
#'
#' # predictions from the model
#' xtest <- seq(0, 1, length = 100)
#' ytest <- sigfun(xtest)
#' pred <- predict(model, xtest)
#' plot(xtest, ytest, type = "l", lty = 2, main = "Kriging predictions")
#' lines(xtest, pred$Phi.test %*% pred$mu, type = "l", col = "blue")
#' lines(xtest, pred$Phi.test %*% pred$xi.map, type = "l", col = "red")
#' legend("right", c("ytest", "mean", "mode"), lty = c(2,1,1),
#'        col = c("black","blue","red"))
#'
#' @importFrom quadprog solve.QP
#' @export
predict.lineqDGP <- function(object, xtest, ...) {
  model <- augment(object)
  if (!is.matrix(xtest))
    xtest <- matrix(xtest, ncol = model$d)

  # passing some terms from the model
  pred <- list()
  class(pred) <- class(model)
  pred$lb <- model$lb
  pred$ub <- model$ub

  # precomputing some term according to Maatouk et al. [2017]
  n_new <- nrow(xtest)
  Phi.test <- basisCompute.lineqDGP(xtest, model$localParam$m, d = 1, model$constrType)
  switch(model$constrType,
         boundedness = {
           pred$Phi.test <- Phi.test
         }, monotonicity = {
           pred$Phi.test  <- cbind(matrix(1, nrow = n_new), Phi.test)
         }, convexity = {
           pred$Phi.test  <- cbind(matrix(1, nrow = n_new), xtest, Phi.test)
         }, {
           stop('constraint "', model$constrType, '" is not supported')
         }
  )

  # # computing the linear system for the QP solver
  # A <- rbind(model$Phi, model$lineqSys$M)
  # b <- rbind(matrix(model$y, ncol = 1), model$lineqSys$g)

  # computing the conditional mean vector and conditional covariance matrix
  # given the interpolation points
  GammaPhit <- model$Gamma %*% t(model$Phi)
  invPhiGammaPhit <- chol2inv(chol(model$Phi %*% GammaPhit))
  pred$mu <- GammaPhit %*% invPhiGammaPhit %*% model$y
  pred$Sigma.xi <- model$Gamma - GammaPhit %*% invPhiGammaPhit %*% t(GammaPhit)

  # computing the conditional mode vector given the interpolation points and
  # the inequality constraints
  # d <- matrix(0, nrow = nrow(model$Gamma))
  # invGamma <- chol2inv(chol(model$Gamma))
  # pred$xi.map <- solve.QP(invGamma, d, t(A), b, nrow(model$x))$solution
  if (min(eigen(pred$Sigma.xi, symmetric = TRUE)$values) <= 0) # numerical stability
    pred$Sigma.xi <- pred$Sigma.xi + 1e-6*diag(nrow(pred$Sigma.xi))
  invSigma <- chol2inv(chol(pred$Sigma.xi))
  pred$xi.map <- solve.QP(invSigma, t(pred$mu) %*% invSigma,
                          t(model$lineqSys$M), model$lineqSys$g)$solution
  return(pred)
}

#' @title Simulation Method for the \code{"lineqDGP"} S3 Class
#' @description Simulation method for the \code{"lineqDGP"} S3 class.
#' @param object an object with class \code{"lineqDGP"}.
#' @param nsim	the number of simulations.
#' @param seed see \code{\link{simulate}}.
#' @param xtest a vector (or matrix) with the test input design.
#' @param ... further arguments passed to or from other methods.
#' @return An object with the simulations of \code{"lineqDGP"} models.
#' \item{x}{A vector (or matrix) with the training input design.}
#' \item{y}{The training output vector at \code{x}.}
#' \item{xtest}{A vector (or matrix) with the test input design.}
#' \item{Phi.test}{A matrix corresponding to the hat basis functions evaluated at \code{xtest}.
#' The basis functions are indexed by rows.}
#' \item{xi.sim}{Posterior sample-path of the finite-dimensional Gaussian vector.}
#' \item{ysim}{Posterior sample-path of the observed GP. Note: \code{ysim = Phi.test \%*\% xi.sim}.}
#'
#' @details The posterior sample-path of the finite-dimensional GP with inequality constraints is computed
#' according to (Maatouk and Bay, 2017).
#' @seealso \code{\link{create.lineqDGP}}, \code{\link{augment.lineqDGP}}, \code{\link{predict.lineqDGP}}
#' @author A. F. Lopez-Lopera.
#'
#' @references Maatouk, H. and Bay, X. (2017),
#' "Gaussian process emulators for computer experiments with inequality constraints".
#' \emph{Mathematical Geosciences},
#' 49(5):557-582.
#' \href{https://link.springer.com/article/10.1007/s11004-017-9673-2}{[link]}
#'
#' @examples
#' # creating the model
#' sigfun <- function(x) return(1/(1+exp(-7*(x-0.5))))
#' x <- seq(0, 1, length = 5)
#' y <- sigfun(x)
#' model <- create(class = "lineqDGP", x, y, constrType = "monotonicity")
#'
#' # updating and expanding the model
#' model$localParam$m <- 30
#' model$kernParam$par <- c(1, 0.2)
#' model <- augment(model)
#'
#' # sampling from the model
#' xtest <- seq(0, 1, length = 100)
#' ytest <- sigfun(xtest)
#' sim.model <- simulate(model, nsim = 50, seed = 1, xtest = xtest)
#' mu <- apply(sim.model$ysim, 1, mean)
#' qtls <- apply(sim.model$ysim, 1, quantile, probs = c(0.05, 0.95))
#' matplot(xtest, t(qtls), type = "l", lty = 1, col = "gray90",
#'         main = "Constrained Kriging model")
#' polygon(c(xtest, rev(xtest)), c(qtls[2,], rev(qtls[1,])), col = 'gray90', border = NA)
#' lines(xtest, ytest, lty = 2)
#' lines(xtest, mu, type = "l", col = "darkgreen")
#' points(x, y, pch = 20)
#' legend("right", c("ytrain","ytest","mean","confidence"), lty = c(NaN,2,1,NaN),
#'        pch = c(20,NaN,NaN,15), col = c("black","black","darkgreen","gray80"))
#'
#' @importFrom stats simulate
#' @export
simulate.lineqDGP <- function(object, nsim = 1, seed = NULL, xtest, ...) {
  model <- augment(object)
  pred <- predict(model, xtest)

  xi.map <- pred$xi.map
  Sigma.xi <- pred$Sigma.xi
  if (min(eigen(pred$Sigma.xi, symmetric = TRUE)$values) < 0)
    Sigma.xi <- pred$Sigma.xi + model$kernParam$nugget*diag(nrow(pred$Sigma.xi))

  # listing control terms
  control <- as.list(unlist(model$localParam$samplingParam))
  control$mvec <- model$localParam$mvec # for HMC
  control$constrType <- model$constrType # for HMC

  # sampling from the truncated multinormal
  tmvPar <- list(mu = xi.map, Sigma = Sigma.xi, lb = pred$lb, ub = pred$ub)
  class(tmvPar) <- model$localParam$sampler
  set.seed(seed)
  xi.sim <- tmvrnorm(tmvPar, nsim, control)

  # passing some terms to the simulated model
  simModel <- list()
  simModel$x <- model$x
  simModel$y <- model$y
  simModel$xtest <- xtest
  simModel$Phi.test <- pred$Phi.test
  simModel$xi.sim <- xi.sim
  simModel$ysim <- pred$Phi.test %*% xi.sim
  class(simModel) <- class(model)
  return(simModel)
}
