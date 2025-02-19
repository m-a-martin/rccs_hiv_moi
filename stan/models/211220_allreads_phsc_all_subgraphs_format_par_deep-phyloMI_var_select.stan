functions {
  int trunc_binomial_rng(int N, real theta){
    int s;
    s = binomial_rng(N, theta);
    while ( s == 0 ){
      s = binomial_rng(N, theta);
    }
    return s;
  }

  // C/O OF Zhi Ling SSHSPH, NUS
  matrix qr_expand(int N) {
  matrix [N,N] A = diag_matrix(rep_vector(1,N));
  matrix [N,N-1] A_qr;
  for(i in 1:N-1) {
    A[N,i] = -1;
    }
    A[N,N] = 0;
    A_qr = qr_Q(A)[,1:(N-1)];
    return A_qr;
  }

  // C/O OF Zhi Ling SSHSPH, NUS
  /** Sample from a multivariate normal with a given diagonal covariance matrix and sum-to-zero constraint
  *
  * @param zeta_2: diagonal of the covariance matrix
  * @param Q: n*(n-1) matrix from QR decomposition, representing the orthonormal basis
    of the n-1 dimensional subspace orthogonal to the all-ones vector in n-space
  * @param coef_low_std: RV in n-1 dimensional space following standard normal distribution
  * @return: RV in n-dimensional space following N(0, D - D 1 (1' D 1)^{-1} 1'D)
    where the diagonal of D is adjusted by n/(n-1) to compensate for the reduced edge variance due to constraints
  */
  vector adj_cov_cholesky_sample(vector zeta_2, matrix Q, vector coef_low_std) {
    int n = size(zeta_2);
    //real adj_factor = n / (n - 1);
    real adj_factor = n * 1.0 / (n - 1);
    matrix[n,n] D = diag_matrix(adj_factor * zeta_2);
    vector[n] one = rep_vector(1.0, n);
    matrix[n,n] D_constrained = D - D * one * one' * D / trace(D);
    matrix[n-1,n-1] S = Q' * D_constrained * Q;
    matrix[n-1,n-1] L = cholesky_decompose(S);
    vector[n-1] coef_low = L * coef_low_std; // Sample from N(0, S)
    vector[n] coef = Q * coef_low; // Sample from N(0, D_constrained)
    return coef;
  }

} 


data{
  // CHARACTERISTICS OF THE DATA
  int<lower=1> N_ind;   // number of individuals
  int<lower=1> N_obs_max; // maximum number of observations per individual
  array [N_ind] int<lower=1, upper=N_obs_max> N_obs; //number of deep-seq observations per individual
  array [N_ind] int<lower=0, upper=N_obs_max> MI_obs; //number of deep-seq observations per individual with identified multiple infection
   // inputs to handle missingness
  int<lower = 0> N_missing; // number of values in the X_mi matrix with missingness
  array[N_missing, 2] int<lower = 0, upper=N_ind> idx_missing; // indices of mixxing values in X_mi matrix
  vector[N_missing] X_mi_missing_std; // standardization values for missing data, non-missing data presumed to be pre-standardized
  matrix[N_missing, 2] missing_prior; // parameters of priors on missing data
  vector[N_missing] missing_max; // maximum value for each missing value
  vector[N_missing] missing_min; // minimum value for each missing value
  // below are optional depending on model
  matrix[N_ind, 5] X_seq; // seq input data 
  matrix[N_ind, 39] X_mi; // mi input data 
}


transformed data{
  int<lower=1> two_times_N_obs_max = 2*N_obs_max;
  int<lower=0, upper=N_ind> N_ind_with_MI_zero;
  int<lower=0, upper=N_ind> N_ind_with_MI_nnzero;
  int j;
  vector[N_missing] missing_norm_constant; // normalization constants for missing values

  // make index of observations with MI_obs==0, and separately MI_obs>0
  N_ind_with_MI_zero = 0;
  for (i in 1:N_ind){
    N_ind_with_MI_zero += (MI_obs[i]==0) ? 1 : 0; 
  }
  
  N_ind_with_MI_nnzero = N_ind - N_ind_with_MI_zero;
  
  array [N_ind_with_MI_zero] int<lower=1, upper=N_ind> MI_zero_idx;
  j = 1;
  for (i in 1:N_ind){
    if (MI_obs[i]==0){
      MI_zero_idx[j] = i;
      j += 1;
    }
  }
  
  array [N_ind_with_MI_nnzero] int<lower=1, upper=N_ind> MI_nnzero_idx;
  j = 1;
  for (i in 1:N_ind){
    if (MI_obs[i]>0){
      MI_nnzero_idx[j] = i;
      j += 1;
    }
  }

  for ( i in 1:N_missing ){
      // is there a built in function for CDF over a range?
      if (missing_min[i] == 0){ 
        missing_norm_constant[i] =  
          lognormal_lcdf(missing_max[i] | 
              missing_prior[i,1],
              missing_prior[i,2]);
      }else if(missing_min[i] > 0){
         missing_norm_constant[i] =
            log(
              lognormal_cdf(missing_max[i] | 
                missing_prior[i,1],
                missing_prior[i,2]) -
              lognormal_cdf(missing_min[i] | 
                missing_prior[i,1],
                missing_prior[i,2]));

      }
    }
}


parameters{
  real logit_prob_mi_baseline;
  real logit_prob_seq_baseline;
  real logit_prob_mi_fpr;
  real<upper = 2.2> logit_prob_mi_fnr;
  vector[N_ind] logit_prob_seq_ind;
  real<lower = 0> logit_prob_seq_ind_sd;
  vector<lower = 0, upper = 1>[N_missing > 0 ? N_missing : 0] X_mi_missing_raw; // todo : does this work with no mi design matrix
  // below are optional depending on model
  sum_to_zero_vector[2] seq_beta1;
  real seq_beta2;
  sum_to_zero_vector[2] seq_beta3;

  vector[1] mi_beta1;
  vector[5] mi_beta2;
  vector[1] mi_beta3;
  vector[1] mi_beta4;
  vector[1] mi_beta5;
  vector[2] mi_beta6;
  vector[1] mi_beta7;
  vector[1] mi_beta8;
  vector[1] mi_beta9;
  vector[1] mi_beta10;
  vector[2] mi_beta11;
  vector[2] mi_beta12;
  vector[1] mi_beta13;
  vector[1] mi_beta14;
  vector[1] mi_beta15;
  vector[1] mi_beta16;

  real<lower=0> mi_r1_global;
  real<lower=0> mi_r2_global;
  vector<lower=0>[39] mi_r1_local;
  vector<lower=0>[39] mi_r2_local;
}


transformed parameters{
  // declare transformed model parameters
  real prob_mi_fpr;
  real prob_mi_fnr;
  vector[N_ind] prob_mi;
  vector[N_ind] prob_seq_1;
  vector[N_ind] prob_seq_mi;
  vector[N_ind] prob_seq_any;
  vector[N_ind] log1m_prob_seq_1;
  vector[N_missing > 0 ? N_missing : 0] X_mi_missing;
  // below are optional depending on model
  vector[5] logit_prob_seq_coeffs;
  real<lower=0> mi_tau;
  vector<lower=0>[39] mi_lambda;
  vector[39] logit_prob_mi_coeffs;
  mi_tau = mi_r1_global * sqrt(mi_r2_global);
  mi_lambda = mi_r1_local .* sqrt(mi_r2_local);
  vector[2] mi_zeta1;
  vector[6] mi_zeta2;
  vector[2] mi_zeta3;
  vector[2] mi_zeta4;
  vector[2] mi_zeta5;
  vector[3] mi_zeta6;
  vector[2] mi_zeta7;
  vector[2] mi_zeta8;
  vector[2] mi_zeta9;
  vector[2] mi_zeta10;
  vector[3] mi_zeta11;
  vector[3] mi_zeta12;
  vector[2] mi_zeta13;
  vector[2] mi_zeta14;
  vector[2] mi_zeta15;
  vector[2] mi_zeta16;

  
  // construct local design matrix with missing data
  // local construction avoids saving all of these objects in the case of missing data
  {
    mi_zeta1 = adj_cov_cholesky_sample(mi_lambda[1:2]^2, qr_expand(2), mi_beta1);
    mi_zeta2 = adj_cov_cholesky_sample(mi_lambda[3:8]^2, qr_expand(6), mi_beta2);
    mi_zeta3 = adj_cov_cholesky_sample(mi_lambda[9:10]^2, qr_expand(2), mi_beta3);
    mi_zeta4 = adj_cov_cholesky_sample(mi_lambda[11:12]^2, qr_expand(2), mi_beta4);
    mi_zeta5 = adj_cov_cholesky_sample(mi_lambda[13:14]^2, qr_expand(2), mi_beta5);
    mi_zeta6 = adj_cov_cholesky_sample(mi_lambda[15:17]^2, qr_expand(3), mi_beta6);
    mi_zeta7 = adj_cov_cholesky_sample(mi_lambda[18:19]^2, qr_expand(2), mi_beta7);
    mi_zeta8 = adj_cov_cholesky_sample(mi_lambda[20:21]^2, qr_expand(2), mi_beta8);
    mi_zeta9 = adj_cov_cholesky_sample(mi_lambda[22:23]^2, qr_expand(2), mi_beta9);
    mi_zeta10 = adj_cov_cholesky_sample(mi_lambda[24:25]^2, qr_expand(2), mi_beta10);
    mi_zeta11 = adj_cov_cholesky_sample(mi_lambda[26:28]^2, qr_expand(3), mi_beta11);
    mi_zeta12 = adj_cov_cholesky_sample(mi_lambda[29:31]^2, qr_expand(3), mi_beta12);
    mi_zeta13 = adj_cov_cholesky_sample(mi_lambda[32:33]^2, qr_expand(2), mi_beta13);
    mi_zeta14 = adj_cov_cholesky_sample(mi_lambda[34:35]^2, qr_expand(2), mi_beta14);
    mi_zeta15 = adj_cov_cholesky_sample(mi_lambda[36:37]^2, qr_expand(2), mi_beta15);
    mi_zeta16 = adj_cov_cholesky_sample(mi_lambda[38:39]^2, qr_expand(2), mi_beta16);
     // combine coefficients into a single vector
    // zetas of size > 1 are already locally shrunk
    // betas of size 1 are not
    logit_prob_mi_coeffs = [mi_zeta1[1],mi_zeta1[2],mi_zeta2[1],mi_zeta2[2],mi_zeta2[3],mi_zeta2[4],mi_zeta2[5],mi_zeta2[6],mi_zeta3[1],mi_zeta3[2],mi_zeta4[1],mi_zeta4[2],mi_zeta5[1],mi_zeta5[2],mi_zeta6[1],mi_zeta6[2],mi_zeta6[3],mi_zeta7[1],mi_zeta7[2],mi_zeta8[1],mi_zeta8[2],mi_zeta9[1],mi_zeta9[2],mi_zeta10[1],mi_zeta10[2],mi_zeta11[1],mi_zeta11[2],mi_zeta11[3],mi_zeta12[1],mi_zeta12[2],mi_zeta12[3],mi_zeta13[1],mi_zeta13[2],mi_zeta14[1],mi_zeta14[2],mi_zeta15[1],mi_zeta15[2],mi_zeta16[1],mi_zeta16[2]]'*mi_tau;
    matrix[1970, 39] X_mi_imputed = X_mi[,];
    X_mi_missing = missing_min + (missing_max - missing_min) .* X_mi_missing_raw;
    for (i in 1:0){
      X_mi_imputed[idx_missing[i,1], idx_missing[i,2]] = X_mi_missing[i] - X_mi_missing_std[i];
    }
    prob_mi = inv_logit( logit_prob_mi_baseline + X_mi_imputed*logit_prob_mi_coeffs );
  }
  
  // calcualte transformed model parameters
  prob_mi_fpr = inv_logit( logit_prob_mi_fpr );
  prob_mi_fnr = inv_logit( logit_prob_mi_fnr );
   // combine coefficients into a single vector
  logit_prob_seq_coeffs = [seq_beta1[1],seq_beta1[2],seq_beta2,seq_beta3[1],seq_beta3[2]]';
  prob_seq_1 = inv_logit( logit_prob_seq_baseline + logit_prob_seq_ind + X_seq*logit_prob_seq_coeffs );
  log1m_prob_seq_1 = log( 1 - prob_seq_1 );
  prob_seq_mi = ( prob_seq_1 .* inv( 2 - prob_seq_1 ) .* ( 1 - prob_mi_fnr ) ) + 
    ( prob_mi_fpr .* 2 .* ( 1 - prob_seq_1 ) ./ ( 2 - prob_seq_1 ) );
  prob_seq_any = 1 - ( 1 - prob_seq_1 ) .* ( 1 - prob_seq_1 );
}


model{
  //priors
  // sequecing success 
  target += normal_lpdf(logit_prob_seq_baseline | 0, 2 );
  target += normal_lpdf(logit_prob_seq_ind | 0, logit_prob_seq_ind_sd);
  target += cauchy_lpdf(logit_prob_seq_ind_sd | 0, 1);
  // multiple infection status
  target += normal_lpdf(logit_prob_mi_baseline | 0, 3.16 );
  // multiple subgraphs
  target += normal_lpdf(logit_prob_mi_fpr | 0, 1);
  target += normal_lpdf(logit_prob_mi_fnr | 0, 1);
  // below are optional depending on model
  seq_beta1 ~ normal(0,2 * inv(sqrt(1 - inv(2))));
  seq_beta2 ~ normal(0,2);
  seq_beta3 ~ normal(0,2 * inv(sqrt(1 - inv(2))));

  mi_beta1 ~ normal(0,1);
  mi_beta2 ~ normal(0,1);
  mi_beta3 ~ normal(0,1);
  mi_beta4 ~ normal(0,1);
  mi_beta5 ~ normal(0,1);
  mi_beta6 ~ normal(0,1);
  mi_beta7 ~ normal(0,1);
  mi_beta8 ~ normal(0,1);
  mi_beta9 ~ normal(0,1);
  mi_beta10 ~ normal(0,1);
  mi_beta11 ~ normal(0,1);
  mi_beta12 ~ normal(0,1);
  mi_beta13 ~ normal(0,1);
  mi_beta14 ~ normal(0,1);
  mi_beta15 ~ normal(0,1);
  mi_beta16 ~ normal(0,1);
  // Wikipedia: the standard Cauchy distribution is the Student's t-distribution with one degree of freedom, and so it may be constructed by any method that constructs the Student's t-distribution.
  // The location-scale t distribution results from compounding a Gaussian distribution (normal distribution) with mean and unknown variance, with an inverse gamma distribution placed over the variance with parameters a = nu/2 and b = nu*mi_tau^2/2.
  mi_r1_local ~ normal(0.0, 1.0);
  mi_r2_local ~ inv_gamma(0.5*2, 0.5*2);
  mi_r1_global ~ normal(0.0, 1.0);
  mi_r2_global ~ inv_gamma(0.5, 0.5);

  // priors for missing data
  if (N_missing > 0){
    target += lognormal_lpdf  (X_mi_missing | 
      missing_prior[,1], missing_prior[,2]);
    // account for trunction of missing data prior distributions
    // due to min and max limits
    target += -missing_norm_constant;
  }

  // likelihood of observing N_obs and MI_obs in an observed window
  for ( i in 1:N_ind )
  {
    target += log_mix(  prob_mi[i],
                        binomial_lpmf(  MI_obs[i] | N_obs[i], prob_seq_mi[i] ) + 
                          binomial_lpmf( N_obs[i] | N_obs_max, prob_seq_any[i] ) - 
                          log1m_exp( two_times_N_obs_max * log1m_prob_seq_1[i] ),
                        binomial_lpmf( MI_obs[i] | N_obs[i], prob_mi_fpr ) + 
                          binomial_lpmf( N_obs[i] | N_obs_max, prob_seq_1[i] ) -
                          log1m_exp( N_obs_max * log1m_prob_seq_1[i] )
                        );
  }
}

 
generated quantities {
  // posterior predictions
  // number of sequenced windows for each sample
  array[N_ind] int N_counts_ppr;
  array[N_ind] int D_counts_ppr;
  vector[N_ind] ind_log_prob_mi;
  real mi_draw; 
  real tmp;

  for (i in 1:N_ind){
    ind_log_prob_mi[i] = log(prob_mi[i]) + 
                          binomial_lpmf(  MI_obs[i] | N_obs[i], prob_seq_mi[i] ) + 
                          binomial_lpmf( N_obs[i] | N_obs_max, prob_seq_any[i] ) - 
                          log1m_exp( two_times_N_obs_max * log1m_prob_seq_1[i] );
    ind_log_prob_mi[i] -= log_mix(  prob_mi[i],
                        binomial_lpmf(  MI_obs[i] | N_obs[i], prob_seq_mi[i] ) + 
                          binomial_lpmf( N_obs[i] | N_obs_max, prob_seq_any[i] ) - 
                          log1m_exp( two_times_N_obs_max * log1m_prob_seq_1[i] ),
                        binomial_lpmf( MI_obs[i] | N_obs[i], prob_mi_fpr ) + 
                          binomial_lpmf( N_obs[i] | N_obs_max, prob_seq_1[i] ) -
                          log1m_exp( N_obs_max * log1m_prob_seq_1[i] ) );
    tmp = exp(ind_log_prob_mi[i]);
    // avoids edge case error 
    // tmp should never be > 1, but sometimes when ind_log_prob_mi is 
    // very close to 0, stan seems to think it is and bernoulli_rng
    // will throw an error
    if (tmp > 1){
      mi_draw = 1;
    }else{
      mi_draw = bernoulli_rng(tmp);
    }
    if ( mi_draw == 1 ){
      N_counts_ppr[i] = trunc_binomial_rng(N_obs_max, prob_seq_any[i]);
      D_counts_ppr[i] = binomial_rng(N_counts_ppr[i], prob_seq_mi[i]);
    }else{
      N_counts_ppr[i] = trunc_binomial_rng(N_obs_max, prob_seq_1[i]);
      D_counts_ppr[i] = binomial_rng(N_counts_ppr[i], prob_mi_fpr);
    }
  }
}
 


