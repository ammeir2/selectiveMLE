library(dplyr)
library(reshape2)
library(selectiveInference)
library(selectiveMLE)
library(ggplot2)

generate.sqrt.Sigma <- function(p,rho,sigsq=1) {
  distances <- as.matrix(dist(1:p,method="manhattan"))

  Sigma <- matrix(rho,ncol=p,nrow=p)^distances*sigsq

  eigSigma <- eigen(Sigma)
  sqrtSigma <- (eigSigma$vectors)%*%diag(sqrt(eigSigma$values))%*%t(eigSigma$vectors)
  return(list(sqrt=sqrtSigma,Sigma=Sigma))
}

generate.regression.data <- function(n,sqrtSigma,numberNonzero,snr=2,ysig=1) {
  #Generating Explanatory variables
  s <- numberNonzero
  X <- matrix(rnorm(n*p), nrow = n, ncol = p)
  X <- X %*% sqrtSigma

  #Sampling coefficients
  beta <- rep(0,p)
  beta[sample(1:p, s)] <- rexp(s,1)*(1-2*rbinom(s,1,0.5))
  mu <- X %*% beta
  y <- as.vector(mu + rnorm(n, sd=sqrt(ysig)))
  return(list(y=y, X=X, beta=beta, mu=mu))
}

# Parameters ----------------
set.seed(140) #!!
set.seed(143)
#set.seed(135)
# n <- 300
# p <- 200
n <- 200
p <- 400
snr <- 0.2
numberNonzero <- 4
rho <- 0.5

# Generating Data ----------------------
Xsqrtsig <- generate.sqrt.Sigma(p, rho, sigsq = 1)$sqrt
s <- 0
while(s < 2 | s > 40) {
  snrEta <- snr
  regData <- generate.regression.data(n, Xsqrtsig, numberNonzero, snr = snr, ysig = 1)
  X <- regData$X[1:n, ]
  X <- apply(X, 2, function(x) (x - mean(x))/sd(x))
  obsMu <- regData$mu[1:n]
  ysig <- sqrt(var(obsMu) / snrEta)
  y <- rnorm(n, mean = obsMu, sd = ysig)
  y <- y - mean(y)

  ysd <- sd(y)
  y <- y / ysd
  lassoFit <- glmnet::cv.glmnet(X, y, standardize = FALSE, intercept = FALSE)
  lambda <- n * lassoFit$lambda.min
  lassoBeta <- as.vector(coef(lassoFit, s = lambda / n))[-1]
  selected <- lassoBeta != 0
  s <- sum(selected)
  yoracle <- regData$mu[1:n] + rnorm(n, sd = ysig)
  yoracle <- yoracle / ysd
}

trueysig <- ysig/ysd
mu <- regData$mu/ysd
true <- coef(lm(mu ~ regData$X[, selected]))[-1]
whichSelected <- which(selected)
naiveFit <- lm(y ~ X[, selected] - 1)
naiveBeta <- coef(naiveFit)

Xm <- X[, selected]
lassoysig <- sd((y - Xm %*% lassoBeta[lassoBeta != 0]))

# selectiveInference (Polyhedral CIs) ------------------
selectiveFit <- NULL
selectiveysig <- lassoysig * sqrt(n / (n - s - 1))
selectiveFit <- fixedLassoInf(X, y, alpha = 0.05,
                              coef(lassoFit, s = "lambda.min")[-1],
                              lambda = lassoFit$lambda.min * n,
                              sigma = selectiveysig)

# Conditional MLE ----------------
assumeConvergence <- 1400
mle <- lassoMLE(y, X, lambda = "lambda.min",
                 ysig = NULL, lassoFit = lassoFit,
                 delay = 20, assumeConvergence = assumeConvergence,
                 iterations = 2000, stepRate = 0.8)

conditional <- mle$conditionalBeta
# Solution Path -------------------------
solutionPath <- mle$solution_path[1:assumeConvergence, ]
solutionPath <- melt(solutionPath)
names(solutionPath) <- c("iter", "variable", "estimate")
ggplot(solutionPath) +
  geom_line(aes(x = iter, y = estimate, col = factor(variable))) +
  theme_bw() +
  geom_hline(yintercept = 0)


# Truncated Samples -------------------------
samples <- mle$coef_sample
samples <- samples[(assumeConvergence):(nrow(samples)), ]
samples <- melt(samples)
names(samples) <- c("replicate", "variable", "sample")
ggplot(samples) +
  geom_density(aes(x = sample, y = ..density..)) +
  facet_wrap(~ variable, scales = "free") +
  theme_bw() + geom_vline(xintercept = 0, linetype = 2)

# Computing Confidence Intervals ----------
conditionalCI <- mle$wald_CI
selectiveCI <- selectiveFit$ci
naiveSD <- summary(naiveFit)$coefficients[, 2]
naiveCI <- cbind(naiveBeta - 1.96 * naiveSD, naiveBeta + 1.96 * naiveSD)

# CI and estimate example --------------------
lassoBeta <- as.vector(coef(lassoFit, s = lambda / n))[-1]
lassoBeta <- lassoBeta[lassoBeta != 0]
offset <- 0.15
naivedat <- data.frame(estimate = naiveBeta, variable = rank(naiveBeta),
                       lCI = naiveCI[, 1], uCI = naiveCI[, 2],
                       method = "naive", offset = offset)
conddat <- data.frame(estimate = conditional, variable = rank(naiveBeta),
                      lCI = conditionalCI[, 1], uCI = conditionalCI[, 2],
                      method = "conditional", offset = 0)
lassodat <- data.frame(estimate = lassoBeta, variable = rank(naiveBeta),
                      lCI = selectiveCI[, 1], uCI = selectiveCI[, 2],
                      method = "lasso", offset = offset * 2)
truedat <- data.frame(estimate = true, variable = rank(naiveBeta),
                       lCI = NA, uCI = NA,
                       method = "true", offset = offset * 3)
forplot <- rbind(naivedat, lassodat, conddat, truedat)
forplot$method <- factor(as.character(forplot$method), levels = c("conditional", "naive", "true", "lasso"))

forci <- forplot
forci <- subset(forplot, !(method %in% c("laasso", "true")))
forci$method <- as.character(forci$method)
forci$method[forci$method == "lasso"] <- "polyhedral"
forci$method[forci$method == "conditional"] <- "conditional-wald"
forci$method <- factor(forci$method, levels = c("conditional-wald", "naive", "polyhedral"))
slack <- 0.01
cimax <- max(conditionalCI) + slack
cimin <- min(conditionalCI) - slack

# pdf("lassoEstimatesCI.pdf",pagecentre=T, width=8,height=3.5 ,paper = "special")
ggplot(subset(forplot)) +
  geom_point(aes(x = variable + offset, y = estimate, shape = method)) +
  geom_segment(data = forci, aes(x = variable + offset, xend = variable + offset,
                   y = pmax(lCI, cimin), yend = pmin(uCI, cimax), col = method, linetype = method)) +
  theme_bw() + geom_hline(yintercept = 0) + xlab("Variable") +
  ylab("Estimates/CIs")
# dev.off()

