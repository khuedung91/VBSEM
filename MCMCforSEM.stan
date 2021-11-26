// Stan program for fitting a simple structural
// equation model via MCMC.
//
// Last changed: 23 NOV 2021.

data {
  int<lower=0> N; // number of individuals
  int<lower=0> K; // number of items
  vector[K] y[N]; // Y matrix of K items for N individuals
  real mu_lambda; //prior for lambda
  real<lower=0> sig2_lambda; // prior for lambda
  real delta_psi[K]; // prior for psi
  real delta_sig2; // prior for sigma^2
  real<lower=0> sigma_nu; // prior for nu
}

parameters {
  vector[N] eta_norm; // normalized eta for each individual
  vector[K] nu; // int for item k
  vector<lower=0>[K-1] lambday; // loading item k, fixing the first to be 1
  real <lower=0> sigma2; // var of the factor
  vector<lower=0>[K] psidiag; // sd of error
}
transformed parameters{
  vector[N] eta;
  real sigma;
  sigma = sqrt(sigma2);
  eta = sigma*eta_norm; 
}

model{
  vector[K] mu[N];
  matrix[K,K] Sigma;
  
  real cond_sd_lambda[K-1];
  vector[K] lambda;
  
  eta_norm ~ normal(0,1) ;
  lambda[1] = 1;
  lambda[2:K] = lambday;
  target += -1.5*log(sigma2) - 0.5*delta_sig2/(sigma2);
  
  for(k in 1:K){
     target += -1.5*log(psidiag[k]) - 0.5*delta_psi[k]/(psidiag[k]);
     nu[k] ~ normal(0,sigma_nu);    
  }
  
  for(k in 1:(K-1) ){
      cond_sd_lambda[k] = sqrt(sig2_lambda*psidiag[k+1]);
      lambday[k] ~ normal(mu_lambda,cond_sd_lambda[k]);
  }
  
  for(i in 1:N){   
    mu[i] = nu + lambda*eta[i];    
  }
  
  Sigma =  diag_matrix(psidiag);
  
  y ~ multi_normal(mu,Sigma); 
}
