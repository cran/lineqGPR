library("lineqGPR")
library("DiceDesign")
library("plot3D")
library("viridis")

rm(list=ls())
set.seed(7)

#### Synthetic data ####
d <- 20
modatan <- function(x, a) return(atan(a*x))
targetFun <- function(x, d) {
  y <- 0
  a <- (1-(1:d)/d)*5
  for (k in 1:d)
    y <- y + modatan(x[, k], a[k])
  return(y)
}

xdesign <- lhsDesign(1e1*d, d, seed = 7)$design
ydesignNoNoise <- targetFun(xdesign, d)
varNoise <- max(range(ydesignNoNoise))*0.01
ydesign <- ydesignNoNoise

#### Constrained model ####
# creating the model
model <- create(class = "lineqAGP", x = xdesign, y = ydesign,
                constrType = rep("monotonicity", d))
model$localParam$m <- rep(5, d)
for (k in 1:d) 
  model$kernParam[[k]]$par <- c(1, 2)
model$nugget <- 1e-7
model$varnoise <- 1e-2 #varNoise

# training the model
# Note: the covariance parameter estimation takes almost 0.5-1 hour
opttime <- proc.time()
model <- lineqGPOptim(model,
                      x0 = unlist(purrr::map(model$kernParam, "par")),
                      eval_f = "logLik",
                      additive = TRUE,
                      opts = list(algorithm = "NLOPT_LD_MMA",
                                  print_level = 3,
                                  ftol_abs = 1e-3,
                                  maxeval = 1e2,
                                  check_derivatives = FALSE),
                      lb = rep(1e-2, 2*d), ub = rep(c(5, 3), d),
                      estim.varnoise = TRUE,
                      bounds.varnoise = c(1e-2, Inf))
opttime <- proc.time() - opttime

# simulating samples from the model
ntest <- 10
xtest  <- matrix(seq(0, 1, length = ntest), nrow = ntest, ncol = d)
ytest <- targetFun(xtest, d)

nsim <- 1e3
sim.model <- simulate(model, nsim = nsim, seed = 7, xtest = xtest)

colormap <- rev(viridis(1e2))
for (k in seq(5, d, 5)) {
  par(mfrow = c(1,2))
  p <- persp3D(x = unique(xtest[, 1]), y = unique(xtest[, 2]),
               z = outer(modatan(xtest[, 1], (1-1/d)*5),
                         modatan(xtest[, k], (1-k/d)*5), "+"),
               xlab = "x1", ylab = paste("x", k, sep = ""), zlab = "y(x1,..., xd)",
               main = "Target function", phi = 20, theta = -30, col = colormap, colkey = FALSE)

  p <- persp3D(x = unique(xtest[, 1]), y = unique(xtest[, 2]),
               z = outer(rowMeans(sim.model$ysim[[1]]),
                         rowMeans(sim.model$ysim[[k]]), "+"),
               xlab = "x1", ylab = paste("x", k, sep = ""), zlab = "u(x1,..., xd)",
               main = "Predictive mean", phi = 20, theta = -30, col = colormap, colkey = FALSE)
}

