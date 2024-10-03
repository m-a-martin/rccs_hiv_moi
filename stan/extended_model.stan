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
  // SEQUENCING SUCCESS MODEL
  int<lower=0> K_seq; // number of predictors of sequencing success
  matrix[N_ind, K_seq] X_seq; // input data to predict sequencing success
  // MULTIPLE INFECTION MODEL
  // inputs with no missing data
  int<lower=1> K_mi_risk_factors; // number of predictors of multiple infection
  matrix[N_ind, K_mi_risk_factors] X_mi_risk_factors; // input data to predict multiple infection status
  // inputs to handle missingness
  int<lower = 0> N_missing; // number of values in the X_mi matrix with missingness
  array[N_missing, 2] int<lower = 0, upper=N_ind> idx_missing; // indices of mixxing values in X_mi matrix
  vector[N_missing] X_mi_missing_std; // standardization values for missing data, non-missing data presumed to be pre-standardized
  matrix[N_missing, 2] missing_prior; // parameters of priors on missing data
  vector[N_missing] missing_max; // maximum value for each missing value
  vector[N_missing] missing_min; // minimum value for each missing value
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
  real logit_prob_MI;
  real logit_prob_seq_baseline;
  real logit_prob_MI_fpr;
  real<upper = 2.2> logit_prob_MI_fnr;
  vector[K_seq] logit_prob_seq_coeffs;
  real<lower = 0> logit_prob_seq_ind_sd;
  vector[N_ind] logit_prob_seq_ind;
  vector[K_mi_risk_factors] logit_prob_MI_coeffs;
  vector<lower = 0, upper = 1>[N_missing > 0 ? N_missing : 0] X_mi_missing_raw;
}

transformed parameters{
  // declare transformed model parameters
  real prob_MI_fpr;
  real prob_MI_fnr;
  vector[N_ind] prob_MI;
  vector[N_ind] prob_seq_1;
  vector[N_ind] prob_seq_MI;
  vector[N_ind] prob_seq_any;
  vector[N_ind] log1m_prob_seq_1;
  vector[N_missing > 0 ? N_missing : 0] X_mi_missing;
  //// construct local design matrix
  {
  // copy over data
  matrix[N_ind, K_mi_risk_factors] X_mi_risk_factors_imputed = X_mi_risk_factors[,];
  // transform missing data based on bounds
  X_mi_missing = missing_min + (missing_max - missing_min) .* X_mi_missing_raw;
  // replace imputed values and standardize 
  // todo vectorize if possible
  for (i in 1:N_missing){
    X_mi_risk_factors_imputed[idx_missing[i,1], idx_missing[i,2]] = X_mi_missing[i] - X_mi_missing_std[i];
  }
  //X_mi_risk_factors_imputed[idx_missing[,1], idx_missing[,2]] = to_matrix(X_mi_missing - X_mi_missing_std);
  // calculate probability of MI
  prob_MI = inv_logit( logit_prob_MI + X_mi_risk_factors_imputed*logit_prob_MI_coeffs);  
  }
  
  // evaluate transformed model parameters
  prob_MI_fpr = inv_logit( logit_prob_MI_fpr );
  prob_MI_fnr = inv_logit( logit_prob_MI_fnr );
  prob_seq_1 = inv_logit( logit_prob_seq_baseline + X_seq*logit_prob_seq_coeffs  + logit_prob_seq_ind);
  log1m_prob_seq_1 = log( 1 - prob_seq_1 );
  prob_seq_MI = ( prob_seq_1 .* inv( 2 - prob_seq_1 ) .* ( 1 - prob_MI_fnr ) ) + 
    ( prob_MI_fpr .* 2 .* ( 1 - prob_seq_1 ) ./ ( 2 - prob_seq_1 ) );
  prob_seq_any = 1 - ( 1 - prob_seq_1 ) .* ( 1 - prob_seq_1 );
}

model{
  //priors
  target += normal_lpdf(logit_prob_MI | 0, 3.16 );
  target += normal_lpdf(logit_prob_seq_baseline | 0, 2 );
  target += normal_lpdf(logit_prob_seq_coeffs | 0, 2 );
  target += normal_lpdf(logit_prob_MI_fpr | 0, 1);
  target += normal_lpdf(logit_prob_MI_fnr | 0, 1);
  target += normal_lpdf(logit_prob_seq_ind | 0, logit_prob_seq_ind_sd);
  target += cauchy_lpdf(logit_prob_seq_ind_sd | 0, 1);
  target += normal_lpdf(logit_prob_MI_coeffs | 0,  1);
  
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
    target += log_mix(  prob_MI[i],
                        binomial_lpmf(  MI_obs[i] | N_obs[i], prob_seq_MI[i] ) + 
                          binomial_lpmf( N_obs[i] | N_obs_max, prob_seq_any[i] ) - 
                          log1m_exp( two_times_N_obs_max * log1m_prob_seq_1[i] ),
                        binomial_lpmf( MI_obs[i] | N_obs[i], prob_MI_fpr ) + 
                          binomial_lpmf( N_obs[i] | N_obs_max, prob_seq_1[i] ) -
                          log1m_exp( N_obs_max * log1m_prob_seq_1[i] )
                        );
  }
}

 

