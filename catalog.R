# Test & Roll: Profit Maximizing A/B Tests

# Meta-analysis and test & roll design for catalog example

# Elea McDonnell Feit, eleafeit@gmail.com
# 18 July 2019

library(rstan)
library(MASS)
source("nn_functions.R")
options(mc.cores = parallel::detectCores())

# Generate synthetic data -----
generate_syn_expt <- function(nexpt, nobs, nholdout, mu1, Delta, sigma, s, omega) {
  t <- rnorm(nexpt, mu1, omega)
  m1 <- rnorm(nexpt, t, sigma[1])
  m2 <- rnorm(nexpt, t + Delta, sigma[2])
  expt <- rep(1:nexpt, each=nobs) 
  treat <- rep(c(rep(0, nholdout), rep(1, nobs-nholdout)), 30)
  Y <- (1 - treat) * rnorm(nobs*nexpt, mean=rep(m1, each=nobs), sd=s[1]) + 
             treat * rnorm(nobs*nexpt, mean=rep(m2, each=nobs), sd=s[2])
  treat <- treat + 1
  data.frame(Y, expt, treat)
}
set.seed(20190627)
data <- generate_syn_expt(nexpt=30, nobs=1000, nholdout=150, 
                          mu1=19.39, Delta=10.67, sigma=c(20.97, 13.48), 
                          s=c(87.69, 179.36), omega=27.25) # pop parameters from paper

# Descriptives
means <- aggregate(data$Y, list(data$expt, data$treat), mean)
te <- means[means$Group.2==2,3] - means[means$Group.2==1,3]
lift <- te/means[means$Group.2==1,3]
pdf("catalog_descriptives.pdf", height=5, width=5)
hist(lift, col="darkgray", freq=FALSE, main="", xlab="Lift", ylab="Percent")
hist(te, col="darkgray", freq=FALSE, main="", xlab="Treatment Effects", ylab="Percent")
dev.off()

# Meta-analysis of catalog data -----
d <- list(n=nrow(data), nexpt=max(data$expt), y=data$Y, expt=data$expt, treat=data$treat)
m <- stan(file="catalog_model.stan", data=d, chains=4, iter=2000, warmup=1000)

mypars <- c("s", "mu1", "Delta", "sigma", "omega")
traceplot(m, pars=mypars)
print(m, pars=mypars)
print(m, pars=c("m", "t"))
est <- summary(m)$summary[,"mean"]

# Design Test & Roll  -----
# to replicate the exact results in the paper, use these rounded values
mu = c(19.39, 30.06) # est["mu1"], est["mu1"] + est["Delta"]
sigma = c(20.97, 13.48) # est["sigma[1]"], est["sigma[2]"]
s = c(87.69, 179.36) # est["s[1]"], est["s[2]"]

# Plot of estimated response distributions
pts <- -70:150
pdf("catalog_priors.pdf", width=5, height=5)
plot(pts, dnorm(pts, mean=mu[1], sd=sigma[1]), type="l", lty=3, ylim=c(0, 0.03),
     xlab=expression("Mean Sales ("*m[1]*", "*m[2]*")"), ylab=expression("Prior Density"))
lines(pts, dnorm(pts, mean=mu[2], sd=sigma[2]))
legend("topright", legend=c("Treated", "Control"), lty=c(1, 3), cex=0.8, bty="n")
abline(v=mu[1], lty=3)
abline(v=mu[2])
pts <- -110:130
plot(pts, dnorm(pts, mean=mu[2]-mu[1], sd=sqrt(sigma[1]^2 + sigma[2]^2)), type="l", 
     xlab=expression("Difference in Mean Sales ("*m[2]*" - "*m[1]*")"), ylab="Prior Density")
abline(v=mu[2]-mu[1])
dev.off()

# Baseline example
mu[2] <- mu[2] - 0.8 # shift for media
pnorm(0, mean=mu[2]-mu[1], sd=sigma[2]) # probability of negative profit
N = 100000
set.seed(20030601)

# Profit-maximizing test
(n_star <- test_size_nn(N=N, s, mu, sigma))
(eval <- test_eval_nn(n=n_star, N, s, mu, sigma))
1- eval$profit_rand / eval$profit_perfect

# Compare to NHT
d <- mu[1]*0.25  #d <- 0.011 # 20th percentile  # 2*sigma/sqrt(pi) # mean  
(n_nht <- test_size_nht(s, d))
(eval_nht <- test_eval_nn(n=n_nht, N=N, s=s, mu=mu, sigma=sigma))

# Compare to NHT with finite population correction
(n_fpc <- test_size_nht(s=s, d=d, N=N))
(eval_fpc <- test_eval_nn(n=n_fpc, N=N, s=s, mu=mu, sigma=sigma))

# Compare to equal size test
(n_star_equal <- test_size_nn(N=N, s, mu, sigma))
test_eval_nn(n=n_star_equal,  N=N, s=s, mu=mu, sigma=sigma)

# Compare to Thompson sampling
set.seed(20030601)
system.time((out <- profit_nn_sim(n=n_star, N, s, mu, sigma, TS=TRUE, R=10000)))
out$regret
out$profit

# Plots showing optimal sample sizes for different media costs (shifing m2)
mu2 <- 0:70
out <- NULL
for (i in 1:length(mu2)) {
  out <- rbind(out, test_size_nn(N=N, s=s, mu=c(mu[1], mu2[i]), sigma=sigma))
}
out <- cbind(out, mu2, mu2-mu[1])
pdf("catalog_te_sensitivity.pdf", height=5, width=5)
plot(out[,4], out[,1]+out[,2], type="l", ylim=c(0, 4200), lwd=2,
     xlab=expression("Expected Catalog Treatment Effect ("*mu[2]*" - "*mu[1]*")"), ylab="Test Size")
lines(out[,4], out[,1], lwd=2)
polygon(c(out[1,4], out[,4], out[nrow(out),4]), 
        c(0, out[,1], 0), col='gray90')
polygon(c(out[1,4], out[,4], out[nrow(out):1,4]), 
         c(out[1,1] + out[1,2], out[,1], out[nrow(out):1,1] + out[nrow(out):1,2]), col='gray80')
text(40,700, expression(n[2]*" treated"))
text(-10,200, expression(n[1]*" control"))
dev.off()
