lassoMLE <- function(y, X, lambda = "lambda.min",
                     ysig = NULL,
                     lassoFit = NULL,
                     delay = 50,
                     assumeConvergence = 850,
                     iterations = 1650,
                     stepCoef = 0.001, stepRate = 0.85) {
  # Standartizing variables
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
  Xminus <- X[, !selected]
  A0 <- (1 / lambda) * t(Xminus) %*% (diag(n) - projmat)
  u0mat <- t(Xminus) %*% t(hatmat)
  A0y <- as.numeric(A0 %*% y)

  zeroCov <- A0 %*% t(A0) * ysig^2
  oneCov <- ysig^2 * solve(XmX)
  precision <- solve(oneCov)
  condSigma <- ysig^2 / diag(XmX)

  samp <- naiveBeta
  betaSample <- matrix(0, nrow = iterations, ncol = k)
  betaSample[1, ] <- samp
  mu <- as.vector(t(Xm) %*% y)
  estimateMat <- matrix(0, nrow = iterations, ncol = k)
  estimate <- naiveBeta

  singular <- svd(zeroCov)
  sqrtmat <- singular$v %*% (sqrt(singular$d) * t(singular$u))

  zerorank <- min(n - k, p - k)
  zeroMean <- rep(0, p - k)

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
               naive = naiveBeta)

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

  return(result)
}
