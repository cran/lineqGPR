#' @title Hat Basis Functions for \code{"lineqAGP"} Models
#' @description Evaluate the hat basis functions for \code{"lineqAGP"} models.
#' @param x a vector (or matrix) with the input data.
#' @param u a vector (or matrix) with the locations of the knots.
#' @param d a number corresponding to the dimension of the input space.
#' @return A matrix with the hat basis functions. The basis functions are indexed by rows.
#'
#' @section Comments:
#' This function was tested mainly for 1D or 2D input spaces. It could change
#' in future versions for higher dimensions.
#' @author A. F. Lopez-Lopera.
#'
#' @references Lopez-Lopera, A. F., Bachoc, F., Durrande, N., and Roustant, O. (2017),
#' "Finite-dimensional Gaussian approximation with linear inequality constraints".
#' \emph{ArXiv e-prints}
#' \href{https://arxiv.org/abs/1710.07453}{[link]}
#'
#' Maatouk, H. and Bay, X. (2017),
#' "Gaussian process emulators for computer experiments with inequality constraints".
#' \emph{Mathematical Geosciences},
#' 49(5): 557-582.
#' \href{https://link.springer.com/article/10.1007/s11004-017-9673-2}{[link]}
#'
#' @examples
#' x <- seq(0, 1, 1e-3)
#' m <- 5
#' u <- seq(0, 1, 1/(m-1))
#' Phi <- basisCompute.lineqAGP(x, u, d = 1)
#' matplot(Phi, type = "l", lty = 2, main = "Hat basis functions with m = 5")
#'
#' @export
basisCompute.lineqAGP <- function(x, u, d = 1) {
  if (!is.matrix(x) || ncol(x) != d)
    x <- matrix(x, ncol = d)
  if (is.list(u) && length(u) == d) {
    m <- as.integer(lapply(u, length))
  } else if (is.double(u)) {
    u <- matrix(u, ncol = d)
    m <- nrow(u)
  } else {
    m <- nrow(u)
  }

  # precomputing some terms
  n <- nrow(x)
  delta <- 1/(m-1)

  # computing the hat basis functions
  if (d == 1){
    distAbs <- abs(outer(x[, 1]/delta, u[, 1]/delta, "-"))
    idx <- distAbs <= 1
    Phi <- matrix(0, n, m)
    Phi[idx] <- 1 - distAbs[idx]
  } else if (d >= 2){
    PhiList <- list()
    for (k in seq(d)) {
      PhiList[[k]] <- basisCompute.lineqGP(x[, k], u[[k]], d = 1)
    }
  }
  return(Phi)
}

#' @title Linear Systems of Inequalities for \code{"lineqAGP"} Models
#' @description Build the linear system of inequalities for \code{"lineqAGP"} models.
#' @param m the number of linear inequality constraints.
#' @param constrType a character string corresponding to the type of the inequality constraint.
#' Options: "boundedness", "monotonicity", "convexity", "linear"
#' @param l the value (or vector) with the lower bound.
#' @param u the value (or vector) with the upper bound.
#' @param A a matrix containing the structure of the linear equations.
#' @param d the value with the input dimension.
#' @param lineqSysType a character string corresponding to the type of the
#' linear system. Options: \code{twosides}, \code{oneside} (see \code{\link{bounds2lineqSys}} for more details).
#' @param constrIdx for d > 1, a logical vector with the indices of active constrained dimensions.
#' @param rmInf If \code{TRUE}, inactive constraints are removed
#' (e.g. \eqn{-\infty \leq x \leq \infty}{-Inf \le x \le Inf}).
#' @return  A list with the linear system of inequalities: \code{list(A,l,u)} (\code{twosides}) or \code{list(M,g)} (\code{oneside}).
#'
#' @section Comments:
#' This function could change in future versions for more types of inequality
#' constraints in higher dimensions.
#' @seealso \code{\link{bounds2lineqSys}}
#' @author A. F. Lopez-Lopera.
#'
#' @references Lopez-Lopera, A. F., Bachoc, F., Durrande, N., and Roustant, O. (2017),
#' "Finite-dimensional Gaussian approximation with linear inequality constraints".
#' \emph{ArXiv e-prints}
#' \href{https://arxiv.org/abs/1710.07453}{[link]}
#'
#' @examples
#' linSys1 <- lineqAGPSys(m = 5, constrType = "boundedness", l = 0, u = 1, lineqSysType = "twosides")
#' linSys1
#' linSys2 <- lineqAGPSys(m = 5, constrType = "boundedness", l = 0, u = 1, lineqSysType = "oneside")
#' linSys2
#'
#' @export
lineqAGPSys <- function(m = nrow(A),
                        constrType = c("boundedness", "monotonicity",
                                       "convexity", "linear", "none"),
                        l = -Inf, u = Inf,  A = diag(m), d = length(m),
                        lineqSysType = "twosides", constrIdx = seq(length(m)),
                        rmInf = TRUE) {
  constrType <- match.arg(constrType)
  if (constrType != "linear") {
    if (d == 1) {
      switch(constrType,
             boundedness = {
               linSys <- bounds2lineqSys(m, l, u, lineqSysType = lineqSysType, rmInf = rmInf)
             }, monotonicity = {
                 if (length(l) == 1) l <- c(-Inf, rep(l, m-1))
                 if (length(u) == 1) u <- c(Inf, rep(u, m-1))
                 A <- diag(m)
                 diag(A[-1, -ncol(A)]) <- -1
                 linSys <- bounds2lineqSys(nrow(A), l, u, A, lineqSysType, rmInf)
             }, convexity = {
                 if (length(l) == 1) l <- c(-rep(Inf, 2), rep(l, m-2))
                 if (length(u) == 1) u <- c(rep(Inf, 2), rep(u, m-2))
                 A <- diag(m)
                 diag(A[-seq(2), -c(ncol(A)-1,ncol(A))]) <- 1
                 diag(A[-seq(2), -c(1,ncol(A))]) <- -2
                 linSys <- bounds2lineqSys(nrow(A), l, u, A, lineqSysType, rmInf)
             }, none = {
               linSys <- bounds2lineqSys(m, l, u, lineqSysType = lineqSysType, rmInf = rmInf)
             }
      )
    } 
  } else {
    linSys <- bounds2lineqSys(nrow(A), l, u, A, lineqSysType, rmInf)
  }
  return(linSys)
}

#' @title Creation Method for the \code{"lineqAGP"} S3 Class
#' @description Creation method for the \code{"lineqAGP"} S3 class.
#' @param x a vector or matrix with the input data. The dimensions should be indexed by columns.
#' @param y a vector with the output data.
#' @param constrType a character string corresponding to the type of the inequality constraint.
#' Options: "boundedness", "monotonicity", "convexity", "linear";
#' Multiple constraints can be also defined, e.g. \code{constrType = c("boundedness", "monotonicity")}.
#' @return A list with the following elements.
#' \item{x,y,constrType}{see \bold{Arguments}.}
#' \item{d}{a number corresponding to the input dimension.}
#' \item{constrIdx}{for d > 1, a integer vector with the indices of active constrained dimensions.}
#' \item{constrParam}{constraint inequalities for each dimension.}
#' \item{varnoise}{a scalar with noise variance.}
#' \item{localParam}{a list with specific parameters required for \code{"lineqAGP"} models:
#' \code{m} (number of basis functions), \code{sampler}, and \code{samplingParams}.
#' See \code{\link{simulate.lineqAGP}}.}
#' \item{kernParam}{a list with the kernel parameters: \code{par} (kernel parameters), \code{type}, \code{nugget}.
#' See \code{\link{kernCompute}}}
#' \item{bounds}{the limit values if \code{constrType = "boundedness"}.}
#' \item{(Lambda,lb,ub)}{the linear system of inequalities if \code{constrType = "linear"}.}
#'
#' @seealso \code{\link{augment.lineqAGP}}, \code{\link{predict.lineqAGP}}, \code{\link{simulate.lineqAGP}}
#' @author A. F. Lopez-Lopera.
#' @references Lopez-Lopera, A. F., Bachoc, F., Durrande, N., and Roustant, O. (2017),
#' "Finite-dimensional Gaussian approximation with linear inequality constraints".
#' \emph{ArXiv e-prints}
#' \href{https://arxiv.org/abs/1710.07453}{[link]}
#'
#' @examples
#' # creating the model
#' d <- 2
#' fun1 <- function(x) return(4*(x-0.5)^2)
#' fun2 <- function(x) return(2*x)
#' targetFun <- function(x) return(fun1(x[, 1]) + fun1(x[, 2])) 
#' xgrid <- expand.grid(seq(0, 1, 0.01), seq(0, 1, 0.01))
#' ygrid <- targetFun(xgrid)
#' xdesign <- rbind(c(0.5, 0), c(0.5, 0.5), c(0.5, 1), c(0, 0.5), c(1, 0.5))
#' ydesign <- targetFun(xdesign)
#' model <- create(class = "lineqAGP", x = xdesign, y = ydesign,
#'                 constrType = c("convexity", "monotonicity"))
#' str(model)
#' 
#' @method create lineqAGP
#' @export
create.lineqAGP <- function(x, y, constrType) {
  # changing the data as matrices
  if (!is.matrix(x))
    x <- as.matrix(x)
  if (!is.matrix(y) || ncol(y) != 1)
    y <- matrix(y)
  d <- ncol(x) # dimension of the input space
  # creating some lists for the model
  localParam <- list(m = rep(5*length(y), d), sampler = "HMC",
                     samplingParams = c(thinning = 1, burn.in = 1, scale = 0.1))
  constrFlags <- rep(1, d)
  
  kernParam <- vector("list", d)
  constrParam <- vector("list", d)
  for (k in 1:d) {
    kernParam[[k]] <- list(par = c(sigma2 = 1^2, theta = 0.1), type = "gaussian")#, nugget = 1e-7*sd(y))
    switch (constrType[k],
            boundedness = {
              constrParam[[k]]$bounds <- c(lower = min(y) - 0.05*abs(max(y) - min(y)),
                                           upper = max(y) + 0.05*abs(max(y) - max(y)))
            }, monotonicity = {
              constrParam[[k]]$bounds <- c(0, Inf)
            }, convexity = {
              constrParam[[k]]$bounds <- c(0, Inf)
            }, linear = {
              # imposing positiveness constraints by default
              constrParam[[k]]$Lambda <- diag(prod(model$localParam$m))
              constrParam[[k]]$lb <- rep(0, nrow(model$Lambda))
              constrParam[[k]]$ub <- rep(Inf, nrow(model$Lambda))
            }, none = {
              constrParam[[k]]$bounds <- c(-Inf, Inf)
              constrFlags[k] <- 0
            }
    )
  }
  constrIdx <- which(constrFlags == 1)
  # creating the full list for the model
  model <- list(x = x, y = y, constrType = constrType, d = d, nugget = 1e-9,
                constrIdx = constrIdx, constrParam = constrParam,
                varnoise = 0,  localParam = localParam, kernParam = kernParam)
  return(model)
}

#' @title Augmenting Method for the \code{"lineqAGP"} S3 Class
#' @description Augmenting method for the \code{"lineqAGP"} S3 class.
#' @param x an object with class \code{lineqGP}.
#' @param ... further arguments passed to or from other methods.
#' @return An expanded \code{"lineqGP"} object with the following additional elements.
#' \item{Phi}{a matrix corresponding to the hat basis functions.
#' The basis functions are indexed by rows.}
#' \item{Gamma}{the covariance matrix of the Gassian vector \eqn{\boldsymbol{\xi}}{\xi}.}
#' \item{(Lambda,lb,ub)}{the linear system of inequalities.}
#' \item{...}{further parameters passed to or from other methods.}
#'
#' @details Some paramaters of the finite-dimensional GP with linear inequality
#' constraints are computed. Here, \eqn{\boldsymbol{\xi}}{\xi} is a centred Gaussian
#' vector with covariance \eqn{\boldsymbol{\Gamma}}{\Gamma}, s.t.
#' \eqn{\boldsymbol{\Phi} \boldsymbol{\xi} = \boldsymbol{y}}{\Phi \xi = y}
#' (interpolation constraints) and
#' \eqn{\boldsymbol{l} \leq \boldsymbol{\Lambda} \boldsymbol{\xi} \leq \boldsymbol{u}}{lb \le \Lambda \xi \le ub}
#' (inequality constraints).
#' @seealso \code{\link{create.lineqAGP}}, \code{\link{predict.lineqAGP}},
#'          \code{\link{simulate.lineqAGP}}
#' @author A. F. Lopez-Lopera.
#'
#' @references Lopez-Lopera, A. F., Bachoc, F., Durrande, N., and Roustant, O. (2017),
#' "Finite-dimensional Gaussian approximation with linear inequality constraints".
#' \emph{ArXiv e-prints}
#' \href{https://arxiv.org/abs/1710.07453}{[link]}
#' 
#' @examples
#' # creating the model
#' d <- 2
#' fun1 <- function(x) return(4*(x-0.5)^2)
#' fun2 <- function(x) return(2*x)
#' targetFun <- function(x) return(fun1(x[, 1]) + fun1(x[, 2])) 
#' xgrid <- expand.grid(seq(0, 1, 0.01), seq(0, 1, 0.01))
#' ygrid <- targetFun(xgrid)
#' xdesign <- rbind(c(0.5, 0), c(0.5, 0.5), c(0.5, 1), c(0, 0.5), c(1, 0.5))
#' ydesign <- targetFun(xdesign)
#' model <- create(class = "lineqAGP", x = xdesign, y = ydesign,
#'                 constrType = c("convexity", "monotonicity"))
#' 
#' # updating and expanding the model
#' model$localParam$m <- rep(50, d)
#' model$kernParam[[1]]$par <- c(1, 0.2)
#' model$kernParam[[2]]$par <- c(1, 0.2)
#' model$nugget <- 1e-9
#' model$varnoise <- 1e-5
#' model <- augment(model)
#' str(model)
#'
#' @importFrom broom augment
#' @export
augment.lineqAGP<- function(x, ...) {
  model <- x
  if (!("nugget" %in% names(model)))
    model$nugget <- 0
  if ("bounds" %in% names(model)) {
    bounds <- model$bounds
  } else {
    bounds <- c(0, Inf)
  }
  # passing some terms from the model
  x <- model$x
  m <- model$localParam$m

  # computing the kernel matrix for the prior
  u <- Gamma <- vector("list", model$d)
  Phi <- vector("list", model$d)
  for (k in 1:model$d) {
    u[[k]] <- matrix(seq(0, 1, by = 1/(m[k]-1)), ncol = 1) # discretization vector
    Gamma[[k]] <- kernCompute(u[[k]], u[[k]], model$kernParam[[k]]$type,
                              model$kernParam[[k]]$par)
    Phi[[k]] <- basisCompute.lineqGP(x[, k], u[[k]])
  }
  model$u <- u
  model$Gamma <- Gamma
  model$Phi <- Phi

  # precomputing the linear system for the QP solver and MCMC samplers
  M <- g <- vector("list", model$d)
  mvec <- vector("list", model$d)
  for (k in 1:model$d) {
    mvec[[k]] <- c()
    if (model$constrType[k] == "linear") {
      if (!("Lambda" %in% names(model)))
        stop('matrix Lambda is not defined')
      Lambda <- model$constrParam[[k]]$Lambda
      lb <- model$constrParam[[k]]$lb
      ub <- model$constrParam[[k]]$ub
      lsys <- lineqAGPSys(nrow(Lambda), model$constrType[k], lb, ub,
                         Lambda, lineqSysType = "oneside")
      lsys2 <- lineqAGPSys(nrow(Lambda), model$constrType[k], lb, ub,
                          Lambda, rmInf = FALSE)
    } else {
      bounds <- model$constrParam[[k]]$bounds
      lsys <- lineqAGPSys(m[k], model$constrType[k], bounds[1], bounds[2],
                          constrIdx = model$constrIdx, lineqSysType = "oneside")
      lsys2 <- lineqAGPSys(m[k], model$constrType[k], bounds[1], bounds[2],
                           constrIdx = model$constrIdx, rmInf = FALSE)
    }
    # oneside linear structure for QP.solver: M = [Lambda,-Lambda] and g = [-lb,ub]
    M[[k]] <- lsys$M
    g[[k]] <- -matrix(lsys$g)
    # twosides linear structure (Lambda, lb, ub) for MCMC samplers
    model$constrParam[[k]]$Lambda <- lsys2$A
    model$constrParam[[k]]$lb <- lsys2$l
    model$constrParam[[k]]$ub <- lsys2$u
    # extra term required for HMC sampler
    mvec[[k]] <- nrow(lsys2$A)
  }
  # adding the parameters to the model structure
  model$lineqSys$M <- M # for QP solve
  model$lineqSys$g <- g # for QP solve
  model$localParam$mvec <- mvec # for HMC sampler
  return(model)
}

#' @title Prediction Method for the \code{"lineqAGP"} S3 Class
#' @description Prediction method for the \code{"lineqAGP"} S3 class.
#' @param object an object with class \code{"lineqAGP"}.
#' @param xtest a vector (or matrix) with the test input design.
#' @param ... further arguments passed to or from other methods.
#' @return A \code{"lineqAGP"} object with the following elements.
#' \item{Lambda}{a matrix corresponding to the linear set of inequality constraints.}
#' \item{lb}{the lower bound vector of the inequalities constraints.}
#' \item{ub}{the upper bound vector of the inequalities constraints.}
#' \item{Phi.test}{a matrix corresponding to the hat basis functions evaluated
#' at \code{xtest}. The basis functions are indexed by rows.}
#' \item{mu}{the unconstrained GP mean predictor.}
#' \item{Sigma}{the unconstrained GP prediction conditional covariance matrix.}
#' \item{xi.map}{the GP maximum a posteriori (MAP) predictor given the inequality constraints.}
#'
#' @details The posterior paramaters of the finite-dimensional GP with linear inequality
#' constraints are computed. Here, \eqn{\boldsymbol{\xi}}{\xi} is a centred Gaussian
#' vector with covariance \eqn{\boldsymbol{\Gamma}}{\Gamma}, s.t.
#' \eqn{\boldsymbol{\Phi} \boldsymbol{\xi} = \boldsymbol{y}}{\Phi \xi = y}
#' (interpolation constraints) and
#' \eqn{\boldsymbol{l} \leq \boldsymbol{\Lambda} \boldsymbol{\xi} \leq \boldsymbol{u}}{lb \le \Lambda \xi \le ub}
#' (inequality constraints).
#' @seealso \code{\link{create.lineqAGP}}, \code{\link{augment.lineqAGP}},
#'          \code{\link{simulate.lineqAGP}}
#' @author A. F. Lopez-Lopera.
#'
#' @references Lopez-Lopera, A. F., Bachoc, F., Durrande, N., and Roustant, O. (2017),
#' "Finite-dimensional Gaussian approximation with linear inequality constraints".
#' \emph{ArXiv e-prints}
#' \href{https://arxiv.org/abs/1710.07453}{[link]}
#'
#' @examples
#' library(plot3D)
#' # creating the model
#' d <- 2
#' fun1 <- function(x) return(4*(x-0.5)^2)
#' fun2 <- function(x) return(2*x)
#' targetFun <- function(x) return(fun1(x[, 1]) + fun2(x[, 2])) 
#' xgrid <- expand.grid(seq(0, 1, 0.01), seq(0, 1, 0.01))
#' ygrid <- targetFun(xgrid)
#' xdesign <- rbind(c(0.5, 0), c(0.5, 0.5), c(0.5, 1), c(0, 0.5), c(1, 0.5))
#' ydesign <- targetFun(xdesign)
#' model <- create(class = "lineqAGP", x = xdesign, y = ydesign,
#'                 constrType = c("convexity", "monotonicity"))
#' 
#' # updating and expanding the model
#' model$localParam$m <- rep(10, d)
#' model$kernParam[[1]]$type <- "matern52"
#' model$kernParam[[2]]$type <- "matern52"
#' model$kernParam[[1]]$par <- c(1, 0.2)
#' model$kernParam[[2]]$par <- c(1, 0.3)
#' model$nugget <- 1e-9
#' model$varnoise <- 1e-5
#' model <- augment(model)
#'
#' # predictions from the model
#' ntest <- 25
#' xtest  <- cbind(seq(0, 1, length = ntest), seq(0, 1, length = ntest))
#' ytest <- targetFun(xtest)
#' pred <- predict(model, xtest)
#' persp3D(x = unique(xtest[, 1]), y = unique(xtest[, 2]),
#'         z = outer(c(pred$Phi.test[[1]] %*% pred$xi.map[, 1]),
#'                   c(pred$Phi.test[[2]] %*% pred$xi.map[, 2]), "+"),
#'         xlab = "x1", ylab = "x2", zlab = "mode(x1,x2)", zlim = c(0, 3),
#'         phi = 20, theta = -30, alpha = 1, colkey = FALSE)
#' points3D(x = xdesign[,1], y = xdesign[,2], z = ydesign, col = "black", pch = 19, add = TRUE)
#'
#' @importFrom quadprog solve.QP
#' @import plot3D
#' @export
predict.lineqAGP <- function(object, xtest, ...) {
  model <- augment(object)
  if (!is.matrix(xtest))
    xtest <- matrix(xtest, ncol = model$d)

  # passing some terms from the model
  pred <- list()
  class(pred) <- class(model)
  pred$constrParam <- model$constrParam

  # precomputing some terms
  invGamma <- Phi.test <- vector("list", model$d)
  for (k in 1:model$d) {
    Phi.test[[k]] <- basisCompute.lineqGP(xtest[, k], model$u[[k]])
    invGamma[[k]] <- chol2inv(chol(model$Gamma[[k]]))
  }
  pred$Phi.test <- Phi.test
  
  # # computing the conditional mean vector and conditional covariance matrix
  # # given the interpolation points
  hfunBigPhi <- parse(text = paste("cbind(",
                                   paste("model$Phi[[", 1:model$d, "]]", sep = "", collapse = ","),
                                   ")", sep = ""))
  hfunBigInvGamma <- parse(text = paste("bdiag(",
                                        paste("invGamma[[", 1:model$d, "]]", sep = "", collapse = ","),
                                        ")", sep = ""))  
  bigPhi <- eval(hfunBigPhi)
  bigInvGamma <- eval(hfunBigInvGamma)
  invGammaPhitPhi <- chol2inv(chol(model$varnoise*bigInvGamma + t(bigPhi)%*%bigPhi))
  invPhiGammaPhitFull <- (diag(length(model$y)) - bigPhi%*%invGammaPhitPhi%*%t(bigPhi))/model$varnoise
  
  mu <- matrix(0, nrow = model$localParam$m, ncol = model$d)
  Sigma <- vector("list", model$d) 
  for (k in 1:model$d) {
    GammaPhitList <- model$Gamma[[k]] %*% t(model$Phi[[k]])
    temp <- GammaPhitList %*% invPhiGammaPhitFull
    mu[, k] <- temp %*% model$y
    Sigma[[k]] <- model$Gamma[[k]] - temp %*%t(GammaPhitList)
  }
  pred$mu <- mu
  pred$Sigma <- Sigma
  
  xi.map <- matrix(0, nrow = model$localParam$m, ncol = model$d)
  for (k in 1:model$d) {
    # computing the conditional mode vector given the interpolation points
    # and the inequality constraints
    if (min(eigen(pred$Sigma[[k]], symmetric = TRUE)$values) <= 0) # numerical stability
      pred$Sigma[[k]] <- pred$Sigma[[k]] + model$nugget*diag(nrow(pred$Sigma[[k]]))
    invSigma <- chol2inv(chol(pred$Sigma[[k]]))
    xi.map[, k] <- solve.QP(invSigma, t(pred$mu[, k]) %*% invSigma,
                            t(model$lineqSys$M[[k]]), model$lineqSys$g[[k]])$solution
  }
  pred$xi.map <- xi.map
  return(pred)
}

#' @title Simulation Method for the \code{"lineqAGP"} S3 Class
#' @description Simulation method for the \code{"lineqAGP"} S3 class.
#' @param object an object with class \code{"lineqAGP"}.
#' @param nsim	the number of simulations.
#' @param seed see \code{\link{simulate}}.
#' @param xtest a vector (or matrix) with the test input design.
#' @param ... further arguments passed to or from other methods.
#' @return A \code{"lineqAGP"} object with the following elements.
#' \item{x}{a vector (or matrix) with the training input design.}
#' \item{y}{the training output vector at \code{x}.}
#' \item{xtest}{a vector (or matrix) with the test input design.}
#' \item{Phi.test}{a matrix corresponding to the hat basis functions evaluated
#' at \code{xtest}. The basis functions are indexed by rows.}
#' \item{xi.sim}{the posterior sample-path of the finite-dimensional Gaussian vector.}
#' \item{ysim}{the posterior sample-path of the observed GP.
#' Note: \code{ysim = Phi.test \%*\% xi.sim}.}
#'
#' @details The posterior sample-path of the finite-dimensional GP with linear inequality
#' constraints are computed. Here, \eqn{\boldsymbol{\xi}}{\xi} is a centred Gaussian
#' vector with covariance \eqn{\boldsymbol{\Gamma}}{\Gamma}, s.t.
#' \eqn{\boldsymbol{\Phi} \boldsymbol{\xi} = \boldsymbol{y}}{\Phi \xi = y}
#' (interpolation constraints) and
#' \eqn{\boldsymbol{l} \leq \boldsymbol{\Lambda} \boldsymbol{\xi} \leq \boldsymbol{u}}{lb \le \Lambda \xi \le ub}
#' (inequality constraints).
#' @seealso \code{\link{create.lineqAGP}}, \code{\link{augment.lineqAGP}},
#'          \code{\link{predict.lineqAGP}}
#' @author A. F. Lopez-Lopera.
#'
#' @references Lopez-Lopera, A. F., Bachoc, F., Durrande, N., and Roustant, O. (2017),
#' "Finite-dimensional Gaussian approximation with linear inequality constraints".
#' \emph{ArXiv e-prints}
#' \href{https://arxiv.org/abs/1710.07453}{[link]}
#'
#' @examples
#' library(plot3D)
#' # creating the model
#' d <- 2
#' fun1 <- function(x) return(4*(x-0.5)^2)
#' fun2 <- function(x) return(2*x)
#' targetFun <- function(x) return(fun1(x[, 1]) + fun2(x[, 2])) 
#' xgrid <- expand.grid(seq(0, 1, 0.01), seq(0, 1, 0.01))
#' ygrid <- targetFun(xgrid)
#' xdesign <- rbind(c(0.5, 0), c(0.5, 0.5), c(0.5, 1), c(0, 0.5), c(1, 0.5))
#' ydesign <- targetFun(xdesign)
#' model <- create(class = "lineqAGP", x = xdesign, y = ydesign,
#'                 constrType = c("convexity", "monotonicity"))
#' 
#' # updating and expanding the model
#' model$localParam$m <- rep(10, d)
#' model$kernParam[[1]]$type <- "matern52"
#' model$kernParam[[2]]$type <- "matern52"
#' model$kernParam[[1]]$par <- c(1, 0.2)
#' model$kernParam[[2]]$par <- c(1, 0.3)
#' model$nugget <- 1e-9
#' model$varnoise <- 1e-5
#' model <- augment(model)
#'
#' # sampling from the model
#' ntest <- 25
#' xtest  <- cbind(seq(0, 1, length = ntest), seq(0, 1, length = ntest))
#' ytest <- targetFun(xtest)
#' sim.model <- simulate(model, nsim = 1e3, seed = 1, xtest = xtest)
#' persp3D(x = unique(xtest[, 1]), y = unique(xtest[, 2]),
#'         z = outer(rowMeans(sim.model$ysim[[1]]),
#'                   rowMeans(sim.model$ysim[[2]]), "+"),
#'         xlab = "x1", ylab = "x2", zlab = "mode(x1,x2)", zlim = c(0, 3),
#'         phi = 20, theta = -30, alpha = 1, colkey = FALSE)
#' points3D(x = xdesign[,1], y = xdesign[,2], z = ydesign, col = "black", pch = 19, add = TRUE)
#'
#' @importFrom stats simulate
#' @import plot3D
#' @export
simulate.lineqAGP <- function(object, nsim = 1, seed = NULL, xtest, ...) {
  predtime <- proc.time()
  model <- augment(object)
  pred <- predict(model, xtest)
  predtime <- proc.time() - predtime
  
  # computing the transformed conditional mode and covariance matrix given
  # the interpolation points and the inequality constraints
  ysim <- xi.sim <- vector("list", model$d)
  simtime <- proc.time()
  for (k in 1:model$d) {
    Lambda <- pred$constrParam[[k]]$Lambda
    eta.map <- as.vector(Lambda %*% pred$xi.map[, k])
    Sigma.eta <- Lambda %*% pred$Sigma[[k]] %*% t(Lambda)
    if (min(eigen(Sigma.eta, symmetric = TRUE)$values) < 0)
      Sigma.eta <- Sigma.eta + 1e-7*diag(nrow(Sigma.eta))
    
    # listing control terms
    control <- as.list(unlist(model$localParam$samplingParam))
    control$mvec <- model$localParam$mvec[[k]] # for HMC
    control$constrType <- model$constrType[k] # for HMC
    
    # sampling from the truncated multinormal
    tmvPar <- list(mu = eta.map, Sigma = Sigma.eta,
                   lb = pred$constrParam[[k]]$lb,
                   ub = pred$constrParam[[k]]$ub)
    class(tmvPar) <- model$localParam$sampler
    set.seed(seed)
    eta <- tmvrnorm(tmvPar, nsim, control)
    xi.sim[[k]] <- qr.solve(Lambda, eta)
    ysim[[k]] <- pred$Phi.test[[k]] %*% xi.sim[[k]]
  }
  simtime <- proc.time() - simtime
  
  # passing some terms to the simulated model
  simModel <- list()
  simModel$x <- model$x
  simModel$y <- model$y
  simModel$xtest <- xtest
  simModel$Phi.test <- pred$Phi.test
  simModel$xi.map <- pred$xi.map
  simModel$xi.sim <- xi.sim
  simModel$ysim <- ysim
  simModel$predtime <- predtime
  simModel$simtime <- simtime
  class(simModel) <- class(model)
  return(simModel)
}

