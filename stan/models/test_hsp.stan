functions {
  int trunc_binomial_rng(int N, real theta){
    int s;
    s = binomial_rng(N, theta);
    while ( s == 0 ){
      s = binomial_rng(N, theta);
    }
    return s;
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

  sum_to_zero_vector[2] mi_beta1;
  sum_to_zero_vector[6] mi_beta2;
  sum_to_zero_vector[2] mi_beta3;
  sum_to_zero_vector[2] mi_beta4;
  sum_to_zero_vector[2] mi_beta5;
  sum_to_zero_vector[3] mi_beta6;
  sum_to_zero_vector[2] mi_beta7;
  sum_to_zero_vector[2] mi_beta8;
  sum_to_zero_vector[2] mi_beta9;
  sum_to_zero_vector[2] mi_beta10;
  sum_to_zero_vector[3] mi_beta11;
  sum_to_zero_vector[3] mi_beta12;
  sum_to_zero_vector[2] mi_beta13;
  sum_to_zero_vector[2] mi_beta14;
  sum_to_zero_vector[2] mi_beta15;
  sum_to_zero_vector[2] mi_beta16;

  real<lower=0> r1_global;
 real<lower=0> r2_global;
  vector<lower=0>[39] r1_local;
  vector<lower=0>[39] r2_local;
}


transformed parameters{
  // declare transformed model parameters
  real prob_mi_fpr;
  real prob_mi_fnr;
  vector[N_ind] logit_prob_mi;
  vector[N_ind] prob_mi;
  vector[N_ind] prob_seq_1;
  vector[N_ind] prob_seq_mi;
  vector[N_ind] prob_seq_any;
  vector[N_ind] log1m_prob_seq_1;
  vector[N_missing > 0 ? N_missing : 0] X_mi_missing;
  // below are optional depending on model
  vector[5] logit_prob_seq_coeffs;
  real<lower=0> tau;
  vector<lower=0>[39] lambda;
  vector[39] logit_prob_mi_coeffs;
  tau = r1_global * sqrt(r2_global);
  lambda = r1_local .* sqrt(r2_local);
  
  // construct local design matrix with missing data
  // local construction avoids saving all of these objects in the case of missing data
  {
   // combine coefficients into a single vector
  logit_prob_mi_coeffs = [mi_beta1[1],mi_beta1[2],mi_beta2[1],mi_beta2[2],mi_beta2[3],mi_beta2[4],mi_beta2[5],mi_beta2[6],mi_beta3[1],mi_beta3[2],mi_beta4[1],mi_beta4[2],mi_beta5[1],mi_beta5[2],mi_beta6[1],mi_beta6[2],mi_beta6[3],mi_beta7[1],mi_beta7[2],mi_beta8[1],mi_beta8[2],mi_beta9[1],mi_beta9[2],mi_beta10[1],mi_beta10[2],mi_beta11[1],mi_beta11[2],mi_beta11[3],mi_beta12[1],mi_beta12[2],mi_beta12[3],mi_beta13[1],mi_beta13[2],mi_beta14[1],mi_beta14[2],mi_beta15[1],mi_beta15[2],mi_beta16[1],mi_beta16[2]]' .* lambda*tau;
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
  target += normal_lpdf(logit_prob_seq_coeffs | 0, 2 );
  target += normal_lpdf(mi_beta1 | 0, 1);
  target += normal_lpdf(mi_beta2 | 0, 1);
  target += normal_lpdf(mi_beta3 | 0, 1);
  target += normal_lpdf(mi_beta4 | 0, 1);
  target += normal_lpdf(mi_beta5 | 0, 1);
  target += normal_lpdf(mi_beta6 | 0, 1);
  target += normal_lpdf(mi_beta7 | 0, 1);
  target += normal_lpdf(mi_beta8 | 0, 1);
  target += normal_lpdf(mi_beta9 | 0, 1);
  target += normal_lpdf(mi_beta10 | 0, 1);
  target += normal_lpdf(mi_beta11 | 0, 1);
  target += normal_lpdf(mi_beta12 | 0, 1);
  target += normal_lpdf(mi_beta13 | 0, 1);
  target += normal_lpdf(mi_beta14 | 0, 1);
  target += normal_lpdf(mi_beta15 | 0, 1);
  target += normal_lpdf(mi_beta16 | 0, 1);
  // Wikipedia: the standard Cauchy distribution is the Student's t-distribution with one degree of freedom, and so it may be constructed by any method that constructs the Student's t-distribution.
  // The location-scale t distribution results from compounding a Gaussian distribution (normal distribution) with mean and unknown variance, with an inverse gamma distribution placed over the variance with parameters a = nu/2 and b = nu*tau^2/2.
  r1_local ~ normal(0.0, 1.0);
  r2_local ~ inv_gamma(0.5*2, 0.5*2);
  r1_global ~ normal(0.0, 1.0);
  r2_global ~ inv_gamma(0.5, 0.5);

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
 


