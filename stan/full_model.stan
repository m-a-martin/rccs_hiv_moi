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
  int<lower=1> N_ind;   // number of individuals
  int<lower=1> N_obs_max; // maximum number of observations per individual
  int<lower=0> K_seq; // number of predictors of sequencing success
  array [N_ind] int<lower=1, upper=N_obs_max> N_obs; //number of deep-seq observations per individual
  array [N_ind] int<lower=0, upper=N_obs_max> MI_obs; //number of deep-seq observations per individual with identified multiple infection
  matrix[N_ind, K_seq] X_seq; // covariates to predict sequencing success
}
transformed data{
  int<lower=1> two_times_N_obs_max = 2*N_obs_max;
  int<lower=0, upper=N_ind> N_ind_with_MI_zero;
  int<lower=0, upper=N_ind> N_ind_with_MI_nnzero;
  int j;
  
  // make index of observations with MI_obs==0, and separately MI_obs>0
  N_ind_with_MI_zero = 0;
  for (i in 1:N_ind){
    N_ind_with_MI_zero += (MI_obs[i]==0) ? 1 : 0; 
  }
  
  N_ind_with_MI_nnzero = N_ind - N_ind_with_MI_zero;
}
parameters{
  real logit_prob_MI;
  real logit_prob_seq_baseline;
  real logit_prob_MI_fpr;
  real<upper = 2.20> logit_prob_MI_fnr;
  real<lower = 0> logit_prob_seq_ind_sd;
  vector[K_seq] logit_prob_seq_coeffs;
  vector[N_ind] logit_prob_seq_ind;
}
transformed parameters{
  real prob_MI;
  real N_nnzero_times_log_prob_MI;
  real prob_MI_fpr;
  real prob_MI_fnr;
  vector[N_ind] prob_seq_1;
  vector[N_ind] prob_seq_MI;
  vector[N_ind] prob_seq_any;
  vector[N_ind] log1m_prob_seq_1;  
  prob_MI = inv_logit( logit_prob_MI );  
  prob_MI_fpr = inv_logit( logit_prob_MI_fpr );
  prob_MI_fnr = inv_logit( logit_prob_MI_fnr );

  N_nnzero_times_log_prob_MI = N_ind_with_MI_nnzero * log( prob_MI );
  prob_seq_1 = inv_logit( logit_prob_seq_baseline + X_seq*logit_prob_seq_coeffs + logit_prob_seq_ind);
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
  
  // likelihood of observing N_obs and MI_obs in an observed window
 for ( i in 1:N_ind )
  {
    target += log_mix(  prob_MI,
                        binomial_lpmf(  MI_obs[i] | N_obs[i], prob_seq_MI[i] ) + 
                          binomial_lpmf( N_obs[i] | N_obs_max, prob_seq_any[i] ) - 
                          log1m_exp( two_times_N_obs_max * log1m_prob_seq_1[i] ),
                        binomial_lpmf( MI_obs[i] | N_obs[i], prob_MI_fpr ) + 
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
    ind_log_prob_mi[i] = log(prob_MI) + 
                          binomial_lpmf(  MI_obs[i] | N_obs[i], prob_seq_MI[i] ) + 
                          binomial_lpmf( N_obs[i] | N_obs_max, prob_seq_any[i] ) - 
                          log1m_exp( two_times_N_obs_max * log1m_prob_seq_1[i] );
    ind_log_prob_mi[i] -= log_mix(  prob_MI,
                        binomial_lpmf(  MI_obs[i] | N_obs[i], prob_seq_MI[i] ) + 
                          binomial_lpmf( N_obs[i] | N_obs_max, prob_seq_any[i] ) - 
                          log1m_exp( two_times_N_obs_max * log1m_prob_seq_1[i] ),
                        binomial_lpmf( MI_obs[i] | N_obs[i], prob_MI_fpr ) + 
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
      D_counts_ppr[i] = binomial_rng(N_counts_ppr[i], prob_seq_MI[i]);
    }else{
      N_counts_ppr[i] = trunc_binomial_rng(N_obs_max, prob_seq_1[i]);
      D_counts_ppr[i] = binomial_rng(N_counts_ppr[i], prob_MI_fpr);
    }
  }
}
