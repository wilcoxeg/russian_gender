data {
  int N;                        //no trials
  int dv[N];                    //binary 0,1
  int<lower=1> P;               //no fixefs
  int<lower=0> J;               //no subjects
  int<lower=1> n_u;             //no subj ranefs (equals P)
  int<lower=0> K;               //no items
  int<lower=1> n_w;             //no item ranefs
  int<lower=1,upper=J> subj[N]; //subject indicator
  int<lower=1,upper=K> item[N]; //item indicator
  row_vector[P] X[N];           //fixef design matrix
  row_vector[n_u] Z_u[N];       //subj ranef design matrix
  row_vector[n_w] Z_w[N];       //item ranef design matrix
}

parameters {
  vector[P] beta;                 //fixef coefs
  vector<lower=0>[n_u] sigma_u;   //subj ranef sd 
  vector<lower=0>[n_w] sigma_w;   //item ranef sd
  cholesky_factor_corr[n_u] L_u;  //cholesky factor of subj ranef corr matrix
  cholesky_factor_corr[n_w] L_w;  //cholesky factor of item ranef corr matrix
  vector[n_u] z_u[J];             //spherical subj ranef
  vector[n_w] z_w[K];             //spherical item ranef
}

transformed parameters {
  vector[n_u] u[J];             //subj ranefs
  vector[n_w] w[K];             //item ranefs
  {
    matrix[n_u,n_u] Sigma_u;    //subj ranef cov matrix
    matrix[n_w,n_w] Sigma_w;    //item ranef cov matrix
    Sigma_u = diag_pre_multiply(sigma_u,L_u);
    Sigma_w = diag_pre_multiply(sigma_w,L_w);
    for(j in 1:J)
      u[j] = Sigma_u * z_u[j];
    for(k in 1:K)
      w[k] = Sigma_w * z_w[k];
  }
}

model {
  //priors
  beta ~ normal(0,1);    
  sigma_u ~ exponential(2);
  sigma_w ~ exponential(2); 
  L_u ~ lkj_corr_cholesky(2.0);
  L_w ~ lkj_corr_cholesky(2.0);
  for (j in 1:J)
    z_u[j] ~ normal(0,1);
  for (k in 1:K)
    z_w[k] ~ normal(0,1);
  //likelihood
  for (i in 1:N)
    dv[i] ~ bernoulli_logit(X[i] * beta + Z_u[i] * u[subj[i]] + Z_w[i] * w[item[i]]);
}

generated quantities {
  real Gram;
  real Gen;
  real Synt;
  real Lex;
  real Gram_x_Synt;
  real Gram_x_Lex;
  real Gram_x_Gen_x_Synt;
  real Gram_x_Gen_x_Lex;
  vector[N] Predict_rt;

  Gram = ((exp(beta[1] + beta[2]*1/2)) / (1+(exp(beta[1] + beta[2]*1/2)))) - ((exp(beta[1] + beta[2]*(-1/2))) / (1+(exp(beta[1] + beta[2]*(-1/2)))));
  
  Gen = ((exp(beta[1] + beta[3]*1/2)) / (1+(exp(beta[1] + beta[3]*1/2)))) - ((exp(beta[1] + beta[3]*(-1/2))) / (1+(exp(beta[1] + beta[3]*(-1/2)))));
  
  Synt = ((exp(beta[1] + beta[4]*2/3)) / (1+(exp(beta[1] + beta[4]*2/3)))) - ((exp(beta[1] + beta[4]*(-2/3))) / (1+(exp(beta[1] + beta[4]*(-2/3)))));
  
  Lex = ((exp(beta[1] + beta[5]*2/3)) / (1+(exp(beta[1] + beta[5]*2/3)))) - ((exp(beta[1] + beta[5]*(-2/3))) / (1+(exp(beta[1] + beta[5]*(-2/3)))));
  
  Gram_x_Synt = ((exp(beta[1] + beta[6]*2/3)) / (1+(exp(beta[1] + beta[6]*2/3)))) - ((exp(beta[1] + beta[6]*(-2/3))) / (1+(exp(beta[1] + beta[6]*(-2/3)))));
  
  Gram_x_Lex = ((exp(beta[1] + beta[7]*2/3)) / (1+(exp(beta[1] + beta[7]*2/3)))) - ((exp(beta[1] + beta[7]*(-2/3))) / (1+(exp(beta[1] + beta[7]*(-2/3)))));
  
  Gram_x_Gen_x_Synt = ((exp(beta[1] + beta[8]*1/3)) / (1+(exp(beta[1] + beta[8]*1/3)))) - ((exp(beta[1] + beta[8]*(-1/3))) / (1+(exp(beta[1] + beta[8]*(-1/3)))));
  
  Gram_x_Gen_x_Lex = ((exp(beta[1] + beta[9]*1/3)) / (1+(exp(beta[1] + beta[9]*1/3)))) - ((exp(beta[1] + beta[9]*(-1/3))) / (1+(exp(beta[1] + beta[9]*(-1/3)))));
  
  for(i in 1:N)
    Predict_rt[i] = bernoulli_logit_rng(X[i] * beta + Z_u[i] * u[subj[i]] + Z_w[i] * w[item[i]]);

}
