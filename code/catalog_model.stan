// Test & Roll: Profit Maximizing A/B Tests

// Stan code for display example meta-analysis

// Elea McDonnell Feit, efeit@gmail.com
// 18 July 2019


data {
  int<lower=1> n; // number of observations
  int<lower=1> nexpt; // number of experiments (each with multiple arms)
  vector[n] y; // observations
  int<lower=1, upper=nexpt> expt[n]; // indictor for which experiment each obs belongs to
  int<lower=1, upper=2> treat[n]; // indictor for whether an arm is "treatment=2 or "control"=1
}
parameters {
  real m[nexpt, 2]; 
  real t[nexpt]; 
  real<lower=0> s[2];
  real mu1;
  real Delta;
  real<lower=0> omega; 
  real<lower=0> sigma[2];
}
transformed parameters {
  real mu2[nexpt];
  for (i in 1:nexpt) {
    mu2[i] = t[i] + Delta;
  }
}
model {
  // likelihood
  for (i in 1:nexpt) {
    t[i] ~ normal(mu1, omega);
    m[i, 1] ~ normal(t[i], sigma[1]); 
    m[i, 2] ~ normal(mu2[i], sigma[2]);
  }
  for (i in 1:n) {
    y[i] ~ normal(m[expt[i], treat[i]], s[treat[i]]); 
  }
}

