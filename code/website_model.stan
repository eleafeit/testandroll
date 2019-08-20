// Test & Roll: Profit Maximizing A/B Tests

// Stan code for website example meta-analysis

// Elea McDonnell Feit, efeit@gmail.com
// 18 July 2019

data {
  int<lower=1> nexpt; // number of experiments
  real<lower=0,upper=1> y[nexpt,2]; // observed mean response for each arm
  int nobs[nexpt,2]; // sample size for each arm (1 and 2)
}
parameters {
  real<lower=0, upper=1> m[nexpt,2]; 
  real<lower=0, upper=1> t[nexpt];
  real<lower=0, upper=1> mu;
  real<lower=0> sigma; 
  real<lower=0> omega; 
}
model {
  // priors
  mu ~ normal(0.5, 0.1);
  omega ~ normal(0, 0.1);   
  sigma ~ normal(0, 0.1); 
  // likelihood
  for (i in 1:nexpt) {
    t[i] ~ normal(mu, omega);
  	m[i,1] ~ normal(t[i], sigma);
    m[i,2] ~ normal(t[i], sigma);
    y[i,1] ~ normal(m[i,1], sqrt(m[i,1]*(1-m[i,1])/nobs[i,1])); 
	  y[i,2] ~ normal(m[i,2], sqrt(m[i,2]*(1-m[i,2])/nobs[i,2]));
  }
}

