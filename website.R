# Test & Roll: Profit Maximizing A/B Tests

# Meta-analysis and test & roll design for website example

# Elea McDonnell Feit, eleafeit@gmail.com
# 18 July 2019

library(rstan)
source("nn_functions.R")
options(mc.cores = parallel::detectCores())

# Generate synthetic data -----
generate_syn_expt <- function(nexpt, nobs, mu, sigma, omega) {
  nobs <- matrix(nobs, nrow=nexpt, ncol=2) # equal test sizes
  y <- matrix(NA, nrow=nexpt, ncol=2)
  for (e in 1:nexpt) {
    t <- rnorm(1, mean=mu, sd=omega)
    m <- rnorm(2, mean=t, sd=sigma)
    while (min(t,m) < 0 | max(t,m) > 1) { # reject values outside [0,1] 
      t <- rnorm(1, mean=mu, sd=omega)
      m <- rnorm(2, mean=t, sd=sigma)
    }
    y[e, 1] <- mean(rnorm(nobs[e,1], mean=m[1], sd=sqrt(m[1]*(1-m[1]))))
    y[e, 2] <- mean(rnorm(nobs[e,2], mean=m[2], sd=sqrt(m[2]*(1-m[2]))))
  }
  list(nexpt=nexpt, y=y, nobs=nobs)
}
set.seed(17701)
d <- generate_syn_expt(nexpt=50, nobs=10000, 
                        mu=0.676, sigma=0.030, omega=0.199) # pop estimates from paper


# Meta-analysis of website test data -----
m <- stan(file="website_model.stan", data=d, chains=4, iter=2000, seed=19348)

mypars=c("mu", "sigma", "omega")
traceplot(m, pars=mypars) 
print(m, pars=mypars) 
mu <- summary(m, pars="mu")$summary[1]
sigma <- summary(m, pars="sigma")$summary[1]

# Design Test & Roll -----
# to replicate the exact results in the paper, use these rounded values
mu = 0.68
sigma = 0.03
N = 100000

# Prior on effect size
pdf("website.pdf")
sqrt(2)*sigma/sqrt(pi) # mean effect size
plot_prior_effect_nn(mu, sigma, abs=TRUE)

# Profit-maximizing
(n_star <- test_size_nn(N, s=sqrt(mu*(1-mu)), sigma=sigma))
(eval <- test_eval_nn(n=n_star, N=N, s=sqrt(mu*(1-mu)), mu=mu, sigma=sigma))

# Compare to NHT 
d <- 0.68*0.02 # 2% lift 
(pnorm(d, mean=0, sd = sqrt(2)*sigma)-0.5)*2 # percentile of TE distribution
(n_nht <- test_size_nht(s=sqrt(mu*(1-mu)), d=d)) # to match the profit maximizing
(eval_nht <- test_eval_nn(n=rep(n_nht, 2), N=N, s=sqrt(mu*(1-mu)), mu=mu, sigma=sigma))
eval$profit - eval_nht$profit

# Compare to NHT with Finite Population Correction
(n_fpc <- test_size_nht(s=sqrt(mu*(1-mu)), d=d, N=100000)) 
(eval_fpc <- test_eval_nn(n=rep(n_fpc, 2), N=N, s=sqrt(mu*(1-mu)), mu=mu, sigma=sigma))
eval$profit - eval_fpc$profit

# Profit and error rate as a function of n -----
# Profit
n <- c(1:19, 2:19*10, 2:19*100, 2:19*1000, 2:5*10000)
out <- NULL
for (i in 1:length(n)) {
  out <- rbind(out, test_eval_nn(n=c(n[i], n[i]), N=N, s=sqrt(mu*(1-mu)), mu=mu, sigma=sigma))
}
plot(out$n1, out$profit, type="l", 
     ylim=c(out$profit_rand[1], out$profit_perfect[1]),
     xlab=expression("Test Size (n"[1]*"=n"[2]*")"), ylab="Expected Conversions")
abline(v=n_star)
text(n_star, 0.696*N, "n*=2,284", pos=4)
abline(v=n_fpc, col="red", lty=2)
text(n_fpc, 0.683*N, expression(n[FPC]), col="red", pos=2)  # hard code
abline(v=n_nht, col="red", lty=3)
text(n_nht, 0.683*N, expression(n[HT]), col="red", pos=4)  # hard code

# Error rate
plot(out$n1, out$error_rate, type="l", ylim=c(0, 0.5),
     xlab=expression("Test Size (n"[1]*"=n"[2]*")"), ylab="Error Rate")
abline(v=n_star)
text(n_star, 0.13, "n*=2,284", pos=4)
abline(v=n_fpc, col="red", lty=2)
text(n_fpc, 0.3, expression(n[FPC]), col="red", pos=2)
abline(v=n_nht, col="red", lty=3)
text(n_nht, 0.3, expression(n[HT]), col="red", pos=4)

# Sample size sensitivities -----
cex <- 1
ylim <- c(0,25000)

# Optimal sample sizes for different N 
N <- c(100, (1:10000)*1000)
out <- NULL; fpc <- NULL
for (i in 1:length(N)) {
  out <- c(out, test_size_nn(N[i], s=sqrt(mu*(1-mu)), sigma=sigma)[1])
  fpc <-  c(fpc, test_size_nht(s=sqrt(mu*(1-mu)), d=d, N=N[i])) 
}
plot(N/1000000, out, type="l", ylim=ylim,
     xlab="Population Size (N, millions)", ylab="Test Size", cex.lab=cex)
text(N[2000]/1000000, out[2000], "n*", pos=1, cex=cex)
lines(N/1000000, fpc, col="red", lty=2)
text(N[200]/1000000, fpc[200], expression(n[FPC]), col="red", pos=4, cex=cex)
abline(h=n_nht, col="red", lty=3)
text(N[2000]/1000000, n_nht, expression(n[HT]), col="red", pos=3, cex=cex)

# Optimal sample sizes for different sigma
N=100000
sigma <- seq(0.001, 0.19, 0.001)
d <- qnorm(0.5+0.125, mean=0, sd=sqrt(2)*sigma) #25th percentile of prior on te
out <- NULL; nht <- NULL; fpc <- NULL
for (i in 1:length(sigma)) {
  out <- c(out, test_size_nn(N, s=sqrt(mu*(1-mu)), sigma=sigma[i])[1])
  nht <- c(nht, test_size_nht(s=sqrt(mu*(1-mu)), d=d[i])) 
  fpc <-  c(fpc, test_size_nht(s=sqrt(mu*(1-mu)), d=d[i], N=N)) 
}
plot(sigma, out, type="l", ylim=ylim,
     xlab=expression(paste("StDev of Prior (", sigma, ")", sep="")), 
     ylab="Test Size", cex.lab=1.3, cex.lab=cex)
text(sigma[15], out[15], "n*", pos=1, cex=cex)
lines(sigma, fpc, col="red", lty=2)
text(sigma[30], fpc[30], expression(n[FPC]), col="red", pos=2, cex=cex)
lines(sigma, nht, col="red", lty=3)
text(sigma[30], nht[30], expression(n[HT]), col="red", pos=4, cex=cex)

# Optimal sample sizes for different s
sigma <- 0.03
d <- 0.68*0.02
s <- seq(0.001, 0.7, 0.002)
out <- NULL; nht <- NULL; fpc <- NULL
for (i in 1:length(s)) {
  out <- rbind(out, test_size_nn(N, s=s[i], sigma=sigma))
  nht <- rbind(nht, test_size_nht(s=s[i], d=d))
  fpc <- rbind(fpc, test_size_nht(s=s[i], d=d, N=100000))
}
plot(s, out[,1], type="l", ylim=ylim,
     xlab="StDev of Response (s)", ylab="Test Size", cex.lab=cex)
text(s[250], out[250,1], "n*", pos=3, cex=cex)
lines(s, fpc, type="l", col="red", lty=2)
text(s[300], fpc[300]-700, expression(n[FPC]), col="red", pos=1, cex=cex)
lines(s, nht, type="l", col="red", lty=3)
text(s[250], nht[250], expression(n[HT]), col="red", pos=2, cex=cex)
dev.off()

# Regret distribution and comparison with Thompson Sampling -----
# ***** SLOW RUNNING SIMULATION *****
s <- mu*(1-mu)
pdf("website_ts.pdf", width=5, height=5)
registerDoParallel(cores=detectCores()-1) # use number of coures minus 1
system.time((out <- profit_nn_sim(n=n_star, N=N, s, mu, sigma, TS=TRUE, R=10000)))
out$profit*N
out$regret

cex <- 1
h <- hist(out$regret_draws[,"test_roll"], plot=F, breaks=seq(-0.01,0.07,0.002)) 
h$density = h$counts/sum(h$counts)*100
plot(h, freq=F, xlab="Regret", ylab="Percent", main="Test & Roll",
     xlim=c(-0.01,0.04), ylim=c(0,100), 
     cex=cex, cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
abline(v=mean(out$regret_draws$test_roll), col="red", lty=3)
text(0.005, 80, sprintf("expected\n regret\n %.2f%%", eval$regret*100),
     pos=4,cex=cex, col="red", cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
h = hist(out$regret_draws[,"thom_samp"], plot=F, breaks=seq(-0.01,0.07,0.002)) 
h$density = h$counts/sum(h$counts)*100
plot(h, freq=F, xlab="Regret", ylab="Percent", main="Thompson Sampling",
     xlim=c(-0.01,0.04), ylim=c(0,100),
     cex=cex,cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
abline(v=mean(out$regret_draws$thom_samp), col="red", lty=3)
text(0.005,80, sprintf("expected\n regret\n %.2f%%", mean(out$regret_draws[,"thom_samp"])*100),
     col="red", pos=4, cex=cex, cex.lab=cex, cex.axis=cex, cex.main=cex, cex.sub=cex)
dev.off()


