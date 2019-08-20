# Test & Roll: Profit Maximizing A/B Tests

# Code for sensitivity analysis on regret for Test and Roll versus Thompson Sampling

# Elea McDonnell Feit, eleafeit@gmail.com
# Updated 18 July 2018

rm(list=ls())
setwd("~/repos/profit_max_experiments")
source("nn_functions.R")
set.seed(19980103)
registerDoParallel(cores=detectCores()-1) # use number of coures minus 1

# Baseline values for website example -----
N=100000; mu = 0.63; sigma=0.03; s=sqrt(mu*(1-mu))
R=10000 # draws per condition

# Values to vary -----
Ns <- c(1000000, 100000, 10000, 1000, 100)
ss <- c(0.001, 0.01, 0.1, 1, 10)
sigmas <- c(0.001, 0.01, 0.03, 0.1, 1, 10)
mus <- c(10, 1, 0.63, 0.1, 0)
mudiffs <- c(0, 0.0001, 0.01, 0.1)
Ks <- c(2, 3, 4, 8) 
runs <- length(c(Ns, ss, sigmas, mus, mudiffs, Ks)) # different parameter setings (which each have R draws per condition)

# Setup tables for results -----
exp_regret <- data.frame(matrix(NA, nrow=runs, ncol=13))
colnames(exp_regret) <- c("var", "val", "n1", "n2", "perfect_info", "test_roll", "thom_samp", 
                          "pi_05", "tr_05", "ts_05", "pi_95", "tr_95", "ts_95")
regret <- data.frame(matrix(NA, nrow=runs*R, ncol=4))
colnames(regret) <- c("var", "val", "thom_samp", "test_roll")

# Simulations -----
b <- 0 # indexes rows in the regret table
# Varing K
for (i in seq_along(Ks)) {
  if (Ks[i] == 2) {
    print(system.time((n <- round(test_size_nn(N, s=s, mu=mu, sigma=sigma)))))
  } else {
    #print(system.time((n <- test_size_nn_sim(N=N, s=s, mu=mu, sigma=sigma, K=Ks[i], R=R)$n)))
    if (Ks[i] == 3) n <- rep(1896, 3)
    if (Ks[i] == 4) n <- rep(1882, 4)
    if (Ks[i] == 8) n <- rep(1008, 8)
  }
  print(system.time((out <- profit_nn_sim(n=n, N=N, s=s, mu=mu, sigma=sigma, K=Ks[i], R=R, TS=TRUE))))
  exp_regret[b+1, 1] <- "K"
  exp_regret[b+1, 2:13] <- c(Ks[i], n[1], n[2], as.vector(t(out$regret)))
  regret[1:R + b*R, 1] <- "K"
  regret[1:R + b*R, 2:4] <- cbind(Ks[i], out$regret_draws[,c("thom_samp", "test_roll")])
  b <- b+1
}
# Varying N (population size)
for (i in seq_along(Ns)) {
  n <- ceiling(test_size_nn(Ns[i], s, mu, sigma)) 
  print(system.time((out <- profit_nn_sim(n=n, N=Ns[i], s=s, mu=mu, sigma=sigma, K=2, TS=TRUE, R=R))))
  exp_regret[b+1, 1] <- "N"
  exp_regret[b+1, 2:13] <- c(Ns[i], n[1], n[2], as.vector(t(out$regret)))
  regret[1:R + b*R, 1] <- "N"
  regret[1:R + b*R, 2:4] <- cbind(Ns[i], out$regret_draws[,c("thom_samp", "test_roll")])
  b <- b+1
}
# Varying s
for (i in seq_along(ss)) {
  n <- ceiling(test_size_nn(N, ss[i], mu, sigma)) 
  print(system.time((out <- profit_nn_sim(n, N, ss[i], mu, sigma, K=2, TS=TRUE, R=R))))
  exp_regret[b+1, 1] <- "s"
  exp_regret[b+1, 2:13] <- c(ss[i], n[1], n[2], as.vector(t(out$regret)))
  regret[1:R + b*R, 1] <- "s"
  regret[1:R + b*R, 2:4] <- cbind(ss[i], out$regret_draws[,c("thom_samp", "test_roll")])
  b <- b+1
}
# Varying sigma
for (i in seq_along(sigmas)) {
  n <- ceiling(test_size_nn(N, s, mu, sigmas[i])) 
  print(system.time((out <- profit_nn_sim(n, N, s, mu, sigmas[i], K=2, TS=TRUE, R=R))))
  exp_regret[b+1, 1] <- "sigma"
  exp_regret[b+1, 2:13] <- c(sigmas[i], n[1], n[2], as.vector(t(out$regret)))
  regret[1:R + b*R, 1] <- "sigma"
  regret[1:R + b*R, 2:4] <- cbind(sigmas[i], out$regret_draws[,c("thom_samp", "test_roll")])
  b <- b+1
}
# Varying mu
for (i in seq_along(mus)) {
  n <- ceiling(test_size_nn(N, s, mus[i], sigma)) 
  print(system.time((out <- profit_nn_sim(n, N, s, mus[i], sigma, K=2, TS=TRUE, R=R))))
  exp_regret[b+1, 1] <- "mu"
  exp_regret[b+1, 2:13] <- c(mus[i], n[1], n[2], as.vector(t(out$regret)))
  regret[1:R + b*R, 1] <- "mu"
  regret[1:R + b*R, 2:4] <- cbind(mus[i], out$regret_draws[,c("thom_samp", "test_roll")])
  b <- b+1
}
# Varying mu1 - mu2
for (i in seq_along(mudiffs)) {
  n <- ceiling(test_size_nn(N, s=c(s,s), mu=c(mu+mudiffs[i], mu), sigma=c(sigma, sigma))) 
  print(system.time((out <- profit_nn_sim(n, N, s=c(s,s), mu=c(mu+mudiffs[i], mu), sigma=c(sigma,sigma), K=2, TS=TRUE, R=R))))
  exp_regret[b+1, 1] <- "mudiff"
  exp_regret[b+1, 2:13] <- c(mudiffs[i], n[1], n[2], as.vector(t(out$regret)))
  regret[1:R + b*R, 1] <- "mudiff"
  regret[1:R + b*R, 2:4] <- cbind(mudiffs[i], out$regret_draws[,c("thom_samp", "test_roll")])
  b <- b+1
}

save.image("~/repos/profit_max_experiments/paper_examples/RegretSenR10000.RData")
load("~/repos/profit_max_experiments/paper_examples/RegretSenR10000.RData")

# Summary table -----
library(xtable)
for (i in c("N", "K", "s", "sigma", "mu", "mudiff")) {
  out <- t(exp_regret[exp_regret$var==i, c("val", "test_roll", "thom_samp")])
  colnames(out) <- out[1,]
  out <- out[-1,]
  out <- out*100
  out <- rbind(out, out[1,] - out[2,])
  rownames(out) <- c("Test & Roll", "Thompson Sampling", "Difference")
  print(xtable(out, digits=2))
}

# Plots -----
# Error bar plot (option 1)
regret_plot <- function(data, xlab, ylim=c(-0.03, 0.07), plot_col=c("black", "darkgray")) {
  plot(x=data$val, y=data$test_roll, type="l", col=plot_col[1],
       ylim=ylim, xlab=xlab, ylab="Regret (95% Range)")
  segments(x0=data$val, y0=data$tr_05, y1=data$tr_95, col=plot_col[1])
  lines(x=data$val, y=data$thom_samp, col=plot_col[2])
  segments(x0=data$val*1.05, y0=data$ts_05, y1=data$ts_95, col=plot_col[2])
  legend("topright", lty=1, col=plot_col, bty="n", 
         legend=c("Profit-Max Test & Roll", "Thompson Sampling"))
  abline(h=0)
}

pdf("paper_examples/website_sensitivity_bars.pdf")
options(scipen=999)
regret_plot(exp_regret[exp_regret$var=="N",], "Population Size (N)")
regret_plot(exp_regret[exp_regret$var=="s",], "Response Standard Deviation (s)")
regret_plot(exp_regret[exp_regret$var=="sigma",], "Prior Standard Deviation (sigma)")
regret_plot(exp_regret[exp_regret$var=="mu",], "Mean Response (mu)")
regret_plot(exp_regret[exp_regret$var=="mudiff",], "Difference in Mean Response (mu1 - mu2)")
regret_plot(exp_regret[exp_regret$var=="K",], "Number of Arms (K)")
dev.off()

# Ridgeplots (option 2)
library(ggridges)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)
library(yarrr)
c1 <- "orange"; c2 <- "darkblue"
c1t <- transparent(c1, 0.6); c2t <- transparent(c2, 0.6)
regret_plot <- function(regret, label, l=-0.005, u=0.015) {
  regret %>% filter(!is.na(val)) %>% gather(method, regret, thom_samp:test_roll) %>% mutate(valf=fct_rev(as.factor(val))) %>%
    ggplot(aes(y=valf, x=regret, color=method, fill=method)) + 
      theme(legend.position="none", legend.title = element_blank()) +
      geom_density_ridges(from=l, to=u, scale=1) + 
      scale_color_manual(values = c(c1t, c2t), guide="none") +
      scale_fill_manual(values = c(c1t, c2t), labels = c("Profit-Max Test & Roll", "Thompson Sampling")) +
      scale_discrete_manual("point_color", values = c(c1t, c2t), guide="none") +
      geom_vline(xintercept = 0) +
      guides(fill = guide_legend(override.aes = list(fill = c(c1t, c2t), color = NA, point_color = NA))) +
      labs(y=label, x="Regret")
}  

pdf("paper_examples/website_sensitivity_ridges.pdf", 
    width=3, height=2.5)
options(scipen=999)
regret_plot(regret[regret$var=="N" & regret$val > 100,], "Population Size (N)")
regret_plot(regret[regret$var=="K",], "Number of Arms (K)")
regret_plot(regret[regret$var=="s" & regret$val > 0.001,], "Response St. Dev. (s)")
regret_plot(regret[regret$var=="sigma",], "Prior St. Dev. (sigma)")
regret_plot(regret[regret$var=="mu",], "Prior mean (mu)")
regret_plot(regret[regret$var=="mudiff",], "Prior Difference (mu1 - mu2)")
dev.off()

regret$diff <- regret$test_roll - regret$thom_samp
aggregate(diff ~ var + val, FUN=mean, data=regret)
hist(regret$diff[regret$var=="N" & regret$val==Ns[1]])
