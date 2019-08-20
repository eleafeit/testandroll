# Test & Roll: Profit Maximizing A/B Tests

# Meta-analysis and test & roll design for online display advertising example

# Elea McDonnell Feit, eleafeit@gmail.com
# 18 July 2019

# display_LewisRao2015Retail.csv and functions.R should be in the working directory
source("nn_functions.R")
library(rstan)

# Meta-analysis for first advertiser in Lewis and Rao -----
lr <- read.csv("display_LewisRao2015Retail.csv")
# data taken from tables 1 and 2 of Lewis and Rao (2015)
# note that they do not report the mean for the treatment group or any effects sizes
c <- c(1:3,5:6) # include only advertiser 1 and eliminate exp 4
d1 <- list(nexpt=length(c), nobs=lr$n1[c], ybar=lr$m[c], s=lr$s[c])
m <- stan(file="display_model.stan", data=d1, seed=20030601, 
           iter=10000)
traceplot(m)
summary(m)$summary
# parameters may be estiamted as the mean of the mcmc draws
mu <- summary(m, pars="mu")$summary[1]
sigma <- summary(m, pars="sigma")$summary[1]

# Design test & roll -----
# to replicate the exact results in the paper, use these rounded values for mu and sigma
mu <- 10.36044  
sigma <- 4.39646
margin <- 0.5
(s <- mean(d1$s))  # average of reported standard deviations
(d <- mean(lr$cost[c])/margin) 
N <- 1000000

# Profit-maximizing
# Assume that mu1=mu2, impling that the prior is that 
# the advertising is worth the cost. 
(n <- test_size_nn(N, s, mu, sigma))
(eval <- test_eval_nn(n, N, s, mu, sigma))

# Compare to NHT to detect ROI = 0 versus ROI = -100
(n_nht <- test_size_nht(s=s, d=d))  
test_eval_nn(n_nht, N, s, mu, sigma) # error N > n[1] + n[2]

# Compare to NHT with finite population correction
(n_fpc <- test_size_nht(s=s, d=d, N=N))  
(eval_fpc <- test_eval_nn(c(n_fpc, n_fpc), N, s, mu, sigma))

# Sensitivity: half sigma and N=10,000,000
(n <- test_size_nn(N=100000000, s=s, sigma=sigma/2, mu=mu))
test_eval_nn(n, N, s=s, sigma=sigma/2, mu=mu)

# Other sensitivities not reported in paper -----
# If mu2 > mu1, resulting test is assymetric
(n <- test_size_nn(N, s=c(s,s), mu=c(mu, mu*2), sigma=rep(sigma,2)))
test_eval_nn(n, N, s=c(s,s), mu=c(mu, mu*2), sigma=rep(sigma,2))

# Vary sigma at quantiles of the posterior
draws <- data.frame(extract(m, pars=c("mu", "sigma")))
sigmas <- quantile(draws$sigma, c(0.025, 0.975))
test_size_nn(N, sigma=sigmas[1], s=s)
test_size_nn(N, sigma=sigmas[2], s=s)

