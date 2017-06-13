#' Maximum Likelihood Inference for Models Selected by the LASSO
#'
#' @description \code{lassoMLE} computes the conditional MLE for models
#' selected by the lasso using a stochastic gradient/contrastive divergence
#' approach. Post-selection confidence intervals are computed based on the
#' theoretical approximate distribution of the conditional MLE.
#'
#' @param y response variable.
#'
#' @param X input matrix of dimensions nobs x nvars.
#'
#' @param lambda regularization parameter value to be used with
#' \code{\link[glmnet]{glmnet}}.
#'
#' @param ysig the residual standard error. If \code{NULL} then the
#' cross-validated lasso estimate will be used.
#'
#' @param lassoFit a lasso fit to the data, of the sort returned by
#' \code{\link[glmnet]{glmnet}} or \code{\link[glmnet]{cv.glmnet}}. If
#' NULL, the lasso fit will be computed within the routine.
#'
#' @param delay number of optimization steps to take before decreasing
#' the gradient step size.
#'
#' @param optimSteps number of stochastic gradient steps to take.
#'
#' @param sampSteps number of samples to take in order to compute the confidence
#' intervals.
#'
#' @param stepCoef fixed step size for stochastic gradient.
#'
#' @param stepRate the rate at which to decrease the step size of the
#' stochastic gradient.
#'
#' @param method optimization method, "selected" is faster than "exact"
#' but may not converge to the exact MLE. See description for details.
#'
#' @details The routine computes the conditional MLE for models selected
#' by the lasso as well as post-selection confidence intervals which are
#' based on the approximate distribution of the conditional MLE. The routine
#' uses the \code{\link[glmnet]{cv.glmnet}} to perform model selection and then
#' estimates the regression coefficietns using a stochastic gradient approach.
#'
#' The distinction between the `exact' and `selected' estimation methods is that
#' for exact all the variables relevant to the selection of the lasso model are sampled,
#' including variables related to the set of variables which were not included in the model.
#' Because the number of unselected variables tends to be far larger than the number
#' of variables included in the model, the `selected' method which only samples the refitted
#' regression coefficients tends to be much faster. The two methods tend to produce very similar
#' results.
#'
#'
#'
#' @return \code{lassoMLE} returns an object of class \code{lassoMLE} which
#'   contains the following variables:
#'
#'   * \code{lassoFit} the \code{\link[glmnet]{glmnet}} or
#'   \code{\link[glmnet]{cv.glmnet}} objects describing the solution path
#'   of the lasso.
#'
#'   * \code{activeSet} the models selected by the lasso.
#'
#'   * \code{conditionalBeta} the conditional MLE for the model selected
#'   by the lasso.
#'
#'   * \code{lassoBeta} the lasso coefficients estimates.
#'
#'   * \code{lassoysig} the cross-validated lasso residual standard
#'   error estimate.
#'
#'   * \code{wald_CI} the conditional post-selection confidence intervals.
#'
#'   * \code{solution_path} the solution path of the stochastic gradient method.
#'   This is not the lasso solution path!
#'
#'   * \code{coef_sample} a matrix of samples from the post-selection distribution
#'   of the refitted regression coefficients. These were used to compute the conditional
#'   confidence intervals.
#'
#' @export
#'

lassoMLE <- function(y, X, lambda = "lambda.min",
                      ysig = NULL,
                     lassoFit = NULL,
                     delay = 50,
                     optimSteps = 1000,
                     sampSteps = 2000,
                     stepCoef = 0.001, stepRate = 0.85,
                     method = c("exact", "selected")) {
  assumeConvergence <- optimSteps
  iterations <- optimSteps + sampSteps

  # Checking method ---------------
  if(length(method) > 1) {
    method <- method[1]
  }
  if(!(method %in% c("exact", "selected"))) {
    stop("Method must be either exact or selected!")
  }


  # Standartizing variables ---------------
  lambda <- "lambda.min"
  sdy <- sd(y)
  y <- y/sdy
  meany <- mean(y)
  y <- y - meany
  n <- length(y)
  p <- ncol(X)

  meanX <- colMeans(X)
  sdX <- apply(X, 2, sd)
  for(i in 1:ncol(X)) X[, i] <- (X[, i] - meanX[i]) / sdX[i]

  # Getting the LASSO fit
  if(is.null(lassoFit)) {
    lassoFit <- cv.glmnet(X, y, standardize = FALSE, intercept = FALSE)
  }

  if(lambda == "lambda.min") {
    lambda <- lassoFit$lambda.min * n
  } else if(lambda == "lambda.1se") {
    lambda <- lassoFit$lambda.1se * n
  } else {
    lambda <- lambda
  }

  lassoBeta <- as.vector(coef(lassoFit, s = lambda / n))[-1]
  selected <- lassoBeta != 0
  whichSelected <- which(selected)

  if(all(selected)) {
    method <- "selected"
  }

  lassoysig <- sd((y - X %*% lassoBeta)) * sqrt((n - 1) / (n - sum(selected)))
  naiveysig <- sd(lm(y ~ X[, selected] - 1)$residuals)
  if(is.null(ysig)) {
    ysig <- lassoysig
  }

  ysig <- lassoysig
  k <- sum(selected)
  Xm <- X[, selected]
  Xty <- as.vector(t(Xm) %*% y)
  naiveBeta <- coef(lm(y ~ Xm - 1))
  XmX <- t(Xm) %*% Xm
  XmXinv <- solve(XmX)
  hatmat <- XmXinv %*% t(Xm)
  projmat <- Xm %*% hatmat
  oneCov <- ysig^2 * solve(XmX)
  precision <- solve(oneCov)
  condSigma <- ysig^2 / diag(XmX)

  if(method == "exact") {
    Xminus <- X[, !selected]
    A0 <- (1 / lambda) * t(Xminus) %*% (diag(n) - projmat)
    u0mat <- t(Xminus) %*% t(hatmat)
    A0y <- as.numeric(A0 %*% y)
    zeroCov <- A0 %*% t(A0) * ysig^2
    singular <- svd(zeroCov)
    zerorank <- min(n - k, p - k)
    zeroMean <- rep(0, p - k)
    sqrtmat <- singular$v %*% (sqrt(singular$d) * t(singular$u))
  } else {
    Xminus <- 0
    A0 <- 0
    u0mat <- matrix(0)
    A0y <- as.vector(0)
    zeroCov <- matrix(0)
    singular <- 0
    zerorank <- 0
    zeroMean <- as.vector(0)
    sqrtmat <- matrix(0)
  }

  samp <- naiveBeta
  betaSample <- matrix(0, nrow = iterations, ncol = k)
  betaSample[1, ] <- samp
  mu <- as.vector(t(Xm) %*% y)
  estimateMat <- matrix(0, nrow = iterations, ncol = k)
  estimate <- naiveBeta

  methodExact <- method == "exact"
  lassoBeta <- as.vector(coef(lassoFit, s = lambda / n))[-1][selected]
  lassoSampler(initEst = lassoBeta, initSamp = naiveBeta,  oneCov = oneCov,
               XmX = XmX, XmXinv = XmXinv,
               condSigma = condSigma, lambda = lambda, ysigsq = ysig^2,
               zeroMean = zeroMean, sqrtZero = sqrtmat,
               u0mat = u0mat,
               n, p, nsamp = 1, burnin = 30,
               Xy = Xty,
               estimateMat = estimateMat, sampMat = betaSample,
               delay = delay, stepRate = stepRate, stepCoef = stepCoef,
               gradientBound = 0.02, assumeConvergence = assumeConvergence,
               naive = naiveBeta, methodExact = methodExact)

  conditionalBeta <- estimateMat[nrow(estimateMat), ]

  #conditionalBeta <- pmin(abs(naiveBeta), pmax(0, conditionalBeta * sign(naiveBeta))) * sign(naiveBeta)
  ysig <- sd(Xm %*% conditionalBeta - y)
  ysamp <- betaSample %*% XmX / ysig^2
  forQuantiles <- ysamp[assumeConvergence:iterations, ] / sqrt(n)
  center <- colMeans(forQuantiles)
  variance <- var(forQuantiles)
  coefVar <- solve(variance)
  forQuantiles <- t(apply(forQuantiles, 1, function(x) x - center))
  forQuantiles <- forQuantiles %*% coefVar
  coefDistQuantiles <- t(apply(forQuantiles, 2, function(x) quantile(x, c(0.025, 0.975)))) /sqrt(n)
  trueWald <- cbind(conditionalBeta - coefDistQuantiles[, 2], conditionalBeta - coefDistQuantiles[, 1])
 # roud(cbind(trueWald[, 1], naiveBeta, true, conditionalBeta, trueWald[, 2]), 3)

  # Rescaling parameters and computing intercepts
  lassoBeta <- lassoBeta * sdy
  conditionalBeta <- conditionalBeta * sdy
  lassoBeta <- lassoBeta / sdX[selected]
  conidtionalBeta <- conditionalBeta / sdX[selected]
  ysig <- ysig * sdy
  variance <- variance / sdy^2
  condSD <- sqrt(diag(coefVar))
  conditionalResiduals <- y - as.vector(Xm %*% conditionalBeta)
  conditionalIntercept <- summary(lm(conditionalResiduals ~ 1))$coefficients
  #lassoResiduals <- y - as.vector(Xm %*% lassoBeta)
  lassoIntercept <- summary(lm(conditionalResiduals ~ 1))$coefficients
  #lassoysig <- sd(lassoResiduals)

  result <- list()
  result$lassoFit <- lassoFit
  result$activeSet <- which(selected)
  result$conditionalBeta <- conditionalBeta
  result$lassoBeta <- lassoBeta
  result$lassoysig <- lassoysig
  result$wald_CI <- trueWald
  result$solution_path <- estimateMat
  result$coef_sample <- betaSample

  class(result) <- c("lassoMLE")

  return(result)
}
