/* Fully vectorised Stan code employing Cholesky decomposition for
   McElreath's (2020) multilevel linear model m14.1
   Cf.: Stan User’s Guide and Reference Manual 9.13. Multivariate Priors
     for Hierarchical Models
   Date: Do, 18.06.2020
*/
functions {
}

data {
  /* Dimensions */
  int<lower=1> N;           // num observations
  int<lower=1> M;           // num predictor variables plus one
  int<lower=1> K;           // num groups
  int<lower=1> L;           // num group predictors
    
  /* Observed variables */
  matrix[N, M] X;                 // design matrix
  matrix[K, L] u;                 // matrix of group predictors
  vector[N] wait;                 // outcome variable
  int<lower=1,upper=K> cafe[N];   // group index variable
}

transformed data {
}

parameters {
  /* Unobserved variables */
  matrix[L, M] gamma;             // group coeffs
  real<lower=0> sigma;
  
  /* Cholesky decomposition of covariance matrix */
  vector<lower=0>[M] sigma_cafe;
  cholesky_factor_corr[M] L_Rho;
  matrix[M, K] z_cafe;
}

transformed parameters {
  matrix[K, M] beta;
  //vector[N] eta[K]; *** NOT yet working!!! (cf. p 151 in Manual)***

  /* Correlated varying intercepts and slopes */
  beta = u * gamma + ( diag_pre_multiply( sigma_cafe , L_Rho ) * z_cafe )';
  
  /* Linear model linking expected value of outcome variable to predictors */
  //eta = rows_dot_product( beta[cafe] , X );
}

model {
  /* Fixed priors for unobserved variables */
  target += normal_lpdf( to_vector(gamma) | 0 , 5 );
  target += exponential_lpdf( sigma | 1 );

  /* Fixed priors for Cholesky decomposition of covariance matrix */
  target += exponential_lpdf( sigma_cafe | 1 );
  target += lkj_corr_cholesky_lpdf( L_Rho | 2 );
  target += normal_lpdf( to_vector(z_cafe) | 0 , 1 );

  /* Gauss likelihood for outcome variable */
  target += normal_lpdf( wait | rows_dot_product( beta[cafe] , X ) , sigma );
  //target += normal_lpdf( wait | eta[cafe] , sigma );
}

generated quantities {
  vector[K] a_cafe;
  vector[K] b_cafe;
  matrix[M, M] Rho;
  vector[N] log_lik;

  a_cafe = beta[, 1];
  b_cafe = beta[, 2];
  
  Rho = multiply_lower_tri_self_transpose(L_Rho);
  
  for ( i in 1:N ) {
    log_lik[i] = normal_lpdf( wait[i] | rows_dot_product( beta[cafe] , X )[i] ,
      sigma ); 
  }
}

