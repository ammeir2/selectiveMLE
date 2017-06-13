#' @useDynLib selectiveMLE
#' @importFrom Rcpp sourceCpp

#' @export
truncNormMLE <- function(y, sigma, threshold,
                         cialpha = 0.05,
                         maxiter = 1000,
                         stepRate = 0.7,
                         stepCoef = 1,
                         verbose = TRUE) {
  # Validating input -----------------
  if(ncol(sigma) != nrow(sigma)) {
    stop("sigma must be a symmetric matrix!")
  }

  p <- length(y)
  if(p != ncol(sigma)) {
    stop("length of y must equal the dimension of sigma!")
  }

  if(length(threshold) == 1) {
    threshold <- abs(threshold)
    threshold <- cbind(rep(-threshold, p), rep(threshold, p))
  } else if(length(threshold) == p) {
    threshold <- cbind(-abs(threshold), abs(threshold))
  } else if(all(dim(threshold) == c(p, 2))) {
    threshold <- threshold
  } else {
    stop("threshold must be either a scalar, a vector of size length(y) or
         a matrix of size length(y) X 2")
  }

  selected <- as.integer(as.vector((y < threshold[, 1]) | (y > threshold[, 2])))
  if(all(selected == 0)) {
    stop("at least one coordinate of y must cross the threshold!")
  }

  # Converting to correlations ------------
  sds <- sqrt(diag(sigma))
  sigma <- cov2cor(sigma)
  y <- y / sds
  threshold[, 1] <- threshold[, 1] / sds
  threshold[, 2] <- threshold[, 2] / sds

  # Initializing Values -------------------
  estimate <- y
  samp <- y
  precision <- solve(sigma)
  solutionPath <- matrix(nrow = maxiter, ncol = p)
  solutionPath[1, ] <- estimate
  boolSelected <- selected == 1

  # Optimization Routine --------------------
  if(verbose) print("Starting Optimization")
  for(i in 2:maxiter) {
    sampleMat <- mvtSampler(samp, estimate, selected, threshold,
                            precision, 5, 5, 5)
    samp <- sampleMat[nrow(sampleMat), ]
    expApprox <- colMeans(sampleMat)
    grad <- (precision[boolSelected, , drop = FALSE] %*% (y - expApprox))
    stepSize <- stepCoef / (i - 1)^stepRate
    grad <- grad * stepSize
    grad <- sign(grad) * pmin(abs(grad), 0.01)
    estimate[boolSelected] <- estimate[boolSelected] + grad
    solutionPath[i, ] <- estimate
    if((i %% 100) == 0 & verbose) {
      cat(round(i / maxiter, 2) * 100, "% ")
    }
  }
  if(verbose) cat("\n")

  # Computing Confidence Intervals -----------------
  if(verbose) print("Computing Confidence Intervals")
  mle <- colMeans(solutionPath[max(1, maxiter - 200):maxiter, ])
  sampleMat <- mvtSampler(samp, mle, selected, threshold,
                          precision, max(p * 50, 4000), 500, 25)
  forQuantiles <- apply(sampleMat, 2, function(x) x - mean(x))
  forQuantiles <- forQuantiles %*% precision
  A <- diag(p)
  variance <- var(forQuantiles)
  A[boolSelected, ] <- variance[boolSelected, ]
  Ainv <- solve(A)
  forQuantiles <- forQuantiles %*% Ainv

  quantiles <- apply(forQuantiles[, boolSelected, drop = FALSE], 2,
                     function(x) quantile(x, probs = c(1 - cialpha / 2, cialpha / 2)))
  CI <- matrix(nrow = sum(selected), ncol = 2)
  for(i in 1:nrow(CI)) {
    CI[i, ] <- mle[boolSelected][i] - quantiles[, i]
  }

  # Converting back to original scale ---------------
  mle <- mle[boolSelected] * sds[boolSelected]
  CI[, 1] <- CI[, 1] * sds[boolSelected]
  CI[, 2] <- CI[, 2] * sds[boolSelected]
  for(j in 1:p) {
    solutionPath[, j] <- solutionPath[, j] * sds[j]
    sampleMat[, j] <- sampleMat[, j] * sds[j]
  }

  return(list(mle = mle, CI = CI, solutionPath = solutionPath,
              samples = sampleMat))
}
