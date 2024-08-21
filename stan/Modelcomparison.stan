data {
  int<lower=0> N;           
  vector[N] dv;                
  int<lower=0> J;              
  int<lower=0> K;            
  int<lower=1,upper=J> subj[N]; 
  int<lower=1,upper=K> item[N]; 
  int<lower=1> P;           
  row_vector[P] X[N];
}


parameters {
  vector[P] beta;           
  real<lower=0> sigma_e;       
  real<lower = 0>  tau_u;
  real<lower = 0>  tau_w;
  vector[J] u;
  vector[K] w;
}

model {
  // //priors (vaguely informative)
  target += normal_lpdf(beta[1]| 6, 1);
  target += normal_lpdf(beta[2] | 0, 1);
  target += normal_lpdf(beta[3] | 0, 1);
  target += normal_lpdf(beta[4] | 0, 1);
  target += normal_lpdf(beta[5] | 0, 1);
  target += normal_lpdf(beta[6] | 0, 1);
  target += normal_lpdf(sigma_e | 0, 1)  - normal_lccdf(0 | 0, 1);
  target += normal_lpdf(tau_u | 0, 1)  - normal_lccdf(0 | 0, 1);
  target += normal_lpdf(tau_w | 0, 1)  - normal_lccdf(0 | 0, 1);
  target += normal_lpdf(u | 0, tau_u);
  target += normal_lpdf(w | 0, tau_w);
  // //likelihood
  for (i in 1:N)
    target += lognormal_lpdf(dv[i] | X[i] * beta + u[subj[i]] + w[item[i]], sigma_e);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N)
    log_lik[n] = lognormal_lpdf(dv[n] | X[n] * beta + u[subj[n]] + w[item[n]], sigma_e);
    
}
