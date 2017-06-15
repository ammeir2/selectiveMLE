#' Maximum Likelihood Inference for Selected Normal Means
#'
#' @description \code{truncNormMLE} computes the conditional MLE for subsets of
#' mean vectors of multivariate normal distributions. Conditional estimates are
#' computed via a stochastic gradient/contrastive divergence method and post-selection
#' confidence intervals are computed based on the theoretical apporximate distribution
#' of the conditional MLE.
#'
#' @param y the observed normal vector.
#'
#' @param sigma the covariance matrix of \code{y}.
#'
#' @param threshold a scalar, vector of size \code{length(y)} or matrix of
#' dimensions \code{length(y) X 2}. See description for more details.
#'
#' @param cialpha the confidence intervals will be constructed at a
#' confidence level \code{1 - cialpha}.
#'
#' @param maxiter the number of stochastic gradient steps to take.
#'
#' @param stepRate number of optimization steps to take before decreasing
#' the gradient step size.
#'
#' @param optimSteps number of stochastic gradient steps to take.
#'
#' @param sampSteps the rate at which to decrease the step size of the
#' stochastic gradient.
#'
#' @param stepCoef fixed step size for stochastic gradient.
#'
#' @param verbose whether to report the progess of the optimization routine as
#' it runs.
#'
#'
#' @details The routine computes the conditional MLE for selected normal means.
#' The (coordinate wise) selection rule is assumed to be
#' \code{y[i] > threshold[2, i]  | y[i] < threshold[1, i]}. The function expects
#' the threshold input to be either a scalar, in which case it will be converted into
#' a matrix with values \code{c(-abs(threshold), abs(threshold)} at each row, a vector
#' of size \code{length(y)}, in which case it will be converted to matrix via the formula
#' \code{cbind(-abs(threshold), abs(threshold))}, or a matrix of dimension
#' \code{length(y) X 2}.
#'
#' The default optimization parameter values \code{stepRate} and \code{stepCoef} are
#' quite well optimized
#' and should work well for most problems, as \code{y} is standardized before
#' the optimization starts.
#'
#'
#' @return \code{truncNormMLE} returns an object of class \code{truncNormMLE} which
#'   contains the following variables:
#'
#'   * \code{mle} the conditional MLE.
#'
#'   * \code{CI} the selection adjusted confidence intervals.
#'
#'   * \code{solutionPath} the solution path of the stochastic gradient
#'   method.
#'
#'   * \code{sampleMat} samples from the estimated truncated normal distribution.
#'   These samples are used to compute the post-selection confidence intervals.
#'
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
                            precision, 5, 5, 5, FALSE)
    samp <- sampleMat[nrow(sampleMat), ]
    expApprox <- colMeans(sampleMat)
    grad <- (precision[boolSelected, , drop = FALSE] %*% (y - expApprox))
    stepSize <- stepCoef / (i - 1)^stepRate
    grad <- grad * stepSize
    grad <- sign(grad) * pmin(abs(grad), 0.01)
    estimate[boolSelected] <- estimate[boolSelected] + grad
    solutionPath[i, ] <- estimate
    if((i %% 100) == 0 & verbose) {
      cat(round(i / maxiter, 2) * 100, "% ", sep = "")
    }
  }
  if(verbose) cat("\n")

  # Computing Confidence Intervals -----------------
  if(verbose) print("Computing Confidence Intervals")
  mle <- colMeans(solutionPath[max(1, maxiter - 200):maxiter, ])
  sampleMat <- mvtSampler(samp, mle, selected, threshold,
                          precision, max(p * 50, 4000), 500, 25,
                          verbose)
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

  result <- list(mle = mle, CI = CI, solutionPath = solutionPath,
                 samples = sampleMat)
  class(result) <- "truncNormMLE"

  return(result)
}


