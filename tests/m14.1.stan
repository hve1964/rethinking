data {
  int<lower=1> N;           // num observations
  int<lower=1> M;           // num predictor variables plus one
  int<lower=1> K;           // num groups
  
  int afternoon[N];
  vector[N] wait;
  int<lower=1,upper=K> cafe[N];
}

parameters {
  real a;
  real b;
  real<lower=0> sigma;
  
  vector[K] a_cafe;
  vector[K] b_cafe;
  vector<lower=0>[M] sigma_cafe;
  corr_matrix[M] Rho;
}

model {
  vector[N] mu;
  
  a ~ normal( 5 , 2 );
  b ~ normal( -1 , 0.5 );
  sigma ~ exponential( 1 );

  sigma_cafe ~ exponential( 1 );
  Rho ~ lkj_corr( 2 );

  {
  vector[M] YY[K];
  vector[M] Mu;
  Mu = [ a , b ]';
  for ( j in 1:K ) {
    YY[j] = [ a_cafe[j] , b_cafe[j] ]';
  }
  YY ~ multi_normal( Mu , quad_form_diag(Rho , sigma_cafe) );
  }

  for ( i in 1:N ) {
    mu[i] = a_cafe[cafe[i]] + b_cafe[cafe[i]] * afternoon[i];
  }

  wait ~ normal( mu , sigma );
}

generated quantities {
  vector[N] log_lik;
  vector[N] mu;
  for ( i in 1:N ) {
    mu[i] = a_cafe[cafe[i]] + b_cafe[cafe[i]] * afternoon[i];
  }
  for ( i in 1:N ) {
    log_lik[i] = normal_lpdf( wait[i] | mu[i] , sigma );
  }
}

