library(reshape2)
library(ggplot2)
# Simulation Parameters --------------
set.seed(123)
#set.seed(501)
p <- 100
rho <- 0.3
sigma <- matrix(rho, nrow = p, ncol = p)
diag(sigma) <- 1
threshold <- 1.65
plarge <- 0.2
signalsd <- 2
nullsd <- 0.0000005
cialpha <- 0.05

# Generating Data ------------
ind <- rbinom(p, 1, plarge)
mu <- (1 - ind) * rnorm(p, mean = 0, sd = nullsd) + ind * rnorm(p, mean = 0, sd = signalsd)
mu <- rep(0, p)
mu[1:20] <- rnorm(20, 0, 2)
y <- rep(0, p)
while(all(y < abs(threshold))) {
  y <- as.vector(mvtnorm::rmvnorm(1, mu, sigma))
}
selected <- abs(y) > threshold

# Computing MLE --------------
fit <- truncNormMLE(y, sigma, threshold, cialpha = cialpha,
                    maxiter = 500)

# Plotting Solution Path --------------------
path <- fit$solutionPath
path <- path[, selected]
path <- melt(path)
names(path) <- c("iter", "param", "estimate")
ggplot(path) + geom_line(aes(x = iter, y = estimate, col = factor(param))) +
  geom_hline(yintercept = 0) + theme_bw()

# Plotting Estimates and CIs ------------------------
true <- mu[selected]
offset <- 0.2
naive <- y[selected]
order <- order(naive)
lci <- naive + qnorm(cialpha / 2) * sqrt(diag(sigma)[selected])
uci <- naive + qnorm(1 - cialpha / 2) * sqrt(diag(sigma)[selected])
naivecover <- mean(lci < true & uci > true)
naive <- data.frame(var = order(naive), estimate = naive, lci = lci, uci = uci, method = "naive",
                    offset = offset *  1)
naive$var[order] <- 1:nrow(naive)

conditional <- fit$mle
lci <- fit$CI[, 1]
uci <- fit$CI[, 2]
condcover <- mean(lci < true & uci > true)
conditional <- data.frame(var = naive$var, estimate = conditional, lci = lci, uci = uci, method = "conditional",
                          offset = offset * 0)

true <- data.frame(var = naive$var, estimate = true, lci = NA, uci = NA, method = "true",
                   offset = offset * 2)

forplot <- rbind(naive, conditional, true)
forplot$method <- factor(forplot$method, levels = c("conditional", "naive", "true"))

# pdf("figures/mvtCIexample.pdf",pagecentre=T, width=8,height=4.5 ,paper = "special")
ggplot(forplot, aes(x = var + offset, xend = var + offset)) +
  geom_point(aes(y = estimate, shape = method)) +
  geom_segment(data = subset(forplot, method != "true"), aes(y = lci, yend = uci, col = method, linetype = method)) +
  theme_bw() + geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-threshold, threshold), col = "grey", linetype = 2) +
  xlab("variable #")
# dev.off()

print(c(naivecover, condcover))


