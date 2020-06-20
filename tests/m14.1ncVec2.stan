/* Cf.: T Sorensen, S Hohenstein and S Vasishth (2016)
   Bayesian linear mixed models using Stan: A tutorial for psychologists, 
   linguists, and cognitive scientists
   The Quantitative Methods for Psychology 12 175-200
   Listing 10
   DOI: https://doi.org/10.20982/tqmp.12.3.p175
*/
functions {
}

data {
  /* Dimensions */
  int<lower=1> N;           // num observations
  int<lower=1> M;           // num predictor variables plus one
  int<lower=1> K;           // num groups
    
  /* Observed variables */
  row_vector[M] X[N];             // design matrix
  vector[N] wait;                 // outcome variable
  int<lower=1,upper=K> cafe[N];   // group index variable
}

transformed data {
}

parameters {
  /* Unobserved variables */
  vector[M] beta;               // fixed effects coeffs
  real<lower=0> sigma;
  
  /* Cholesky decomposition of covariance matrix */
  vector<lower=0>[M] sigma_cafe;
  cholesky_factor_corr[M] L_Rho;
  vector[M] z_cafe[K];
}

transformed parameters {
  vector[M] beta_cafe[K];       // varying effects coeffs
  {
    matrix[M, M] Sig_cafe;      // varying effects cov matrix
    Sig_cafe = diag_pre_multiply( sigma_cafe , L_Rho );
    for( j in 1:K ) {
      beta_cafe[j] = Sig_cafe * z_cafe[j];
    }
  }
}

model {
  /* Fixed priors for unobserved variables */
  target += normal_lpdf( beta[1] | 5.0 , 2.0 );
  target += normal_lpdf( beta[2] | -1.0 , 0.5 );
  target += exponential_lpdf( sigma | 1.0 );

  /* Fixed priors for Cholesky decomposition of covariance matrix */
  target += exponential_lpdf( sigma_cafe | 1.0 );
  target += lkj_corr_cholesky_lpdf( L_Rho | 2.0 );
  for ( j in 1:K ) {
    target += normal_lpdf( z_cafe[j] | 0.0 , 1.0 );
  }
 
  /* Gauss likelihood for outcome variable (vectorised) */
  for ( i in 1:N ) {
    target += normal_lpdf( wait[i] | X[i] * (beta + beta_cafe[cafe[i]]) ,
      sigma );
  }
  //target += normal_lpdf( wait | X * (beta + beta_cafe[cafe]) , sigma );
}

generated quantities {
  matrix[M, M] Rho;

  Rho = multiply_lower_tri_self_transpose(L_Rho);
}