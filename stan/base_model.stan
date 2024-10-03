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
}
parameters{
  real logit_prob_MI;
  real logit_prob_seq_baseline;
  vector[K_seq] logit_prob_seq_coeffs;
  vector[N_ind] logit_prob_seq_ind;
  real<lower = 0> logit_prob_seq_ind_sd;
}

transformed parameters{
  real prob_MI;
  real N_nnzero_times_log_prob_MI;
  vector[N_ind] prob_seq_1;
  vector[N_ind] prob_seq_MI;
  vector[N_ind] prob_seq_any;
  vector[N_ind] log1m_prob_seq_1;  
    
  prob_MI = inv_logit( logit_prob_MI );  
  N_nnzero_times_log_prob_MI = N_ind_with_MI_nnzero * log( prob_MI );
  prob_seq_1 = inv_logit( logit_prob_seq_baseline + X_seq*logit_prob_seq_coeffs + logit_prob_seq_ind);
  log1m_prob_seq_1 = log( 1 - prob_seq_1 );
  prob_seq_MI = prob_seq_1 .* inv( 2 - prob_seq_1);
  prob_seq_any = 1 - ( 1 - prob_seq_1 ) .* ( 1 - prob_seq_1 );
}


model{
  //priors
  target += normal_lpdf(logit_prob_MI | 0, 3.16 );
  target += normal_lpdf(logit_prob_seq_baseline | 0, 2 );
  target += normal_lpdf(logit_prob_seq_coeffs | 0, 2 );
  target += normal_lpdf(logit_prob_seq_ind | 0, logit_prob_seq_ind_sd);
  target += cauchy_lpdf(logit_prob_seq_ind_sd | 0, 1);
  
  for ( i in 1:N_ind_with_MI_zero )
  {
    target += log_mix(  prob_MI,
                        binomial_lpmf(  MI_obs[MI_zero_idx[i]] | N_obs[MI_zero_idx[i]], prob_seq_MI[MI_zero_idx[i]] ) + 
                          binomial_lpmf( N_obs[MI_zero_idx[i]] | N_obs_max, prob_seq_any[MI_zero_idx[i]] ) - 
                          log1m_exp( two_times_N_obs_max * log1m_prob_seq_1[MI_zero_idx[i]] ),
                        0 + 
                          binomial_lpmf( N_obs[MI_zero_idx[i]] | N_obs_max, prob_seq_1[MI_zero_idx[i]] ) -
                          log1m_exp( N_obs_max * log1m_prob_seq_1[MI_zero_idx[i]] )
                        );
  }
  for ( i in 1:N_ind_with_MI_nnzero){
    target += log(prob_MI) + 
      binomial_lpmf(MI_obs[MI_nnzero_idx[i]] | N_obs[MI_nnzero_idx[i]], prob_seq_MI[MI_nnzero_idx[i]]) +
      binomial_lpmf( N_obs[MI_nnzero_idx[i]] | N_obs_max, prob_seq_any[MI_nnzero_idx[i]] ) - 
        log1m_exp( two_times_N_obs_max * log1m_prob_seq_1[MI_nnzero_idx[i]] );
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

  for ( i in 1:N_ind_with_MI_zero ){
    ind_log_prob_mi[MI_zero_idx[i]] = log(prob_MI) + 
                        binomial_lpmf(  MI_obs[MI_zero_idx[i]] | N_obs[MI_zero_idx[i]], prob_seq_MI[MI_zero_idx[i]] ) + 
                          binomial_lpmf( N_obs[MI_zero_idx[i]] | N_obs_max, prob_seq_any[MI_zero_idx[i]] ) - 
                          log1m_exp( two_times_N_obs_max * log1m_prob_seq_1[MI_zero_idx[i]] );
    ind_log_prob_mi[MI_zero_idx[i]] -= log_mix(  prob_MI,
                        binomial_lpmf(  MI_obs[MI_zero_idx[i]] | N_obs[MI_zero_idx[i]], prob_seq_MI[MI_zero_idx[i]] ) + 
                          binomial_lpmf( N_obs[MI_zero_idx[i]] | N_obs_max, prob_seq_any[MI_zero_idx[i]] ) - 
                          log1m_exp( two_times_N_obs_max * log1m_prob_seq_1[MI_zero_idx[i]] ),
                        0 + 
                          binomial_lpmf( N_obs[MI_zero_idx[i]] | N_obs_max, prob_seq_1[MI_zero_idx[i]] ) -
                          log1m_exp( N_obs_max * log1m_prob_seq_1[MI_zero_idx[i]] )
                        );
    mi_draw = bernoulli_rng(exp(ind_log_prob_mi[MI_zero_idx[i]]));
    if ( mi_draw == 1 ){
      N_counts_ppr[MI_zero_idx[i]] = trunc_binomial_rng(N_obs_max, prob_seq_any[MI_zero_idx[i]]);
      D_counts_ppr[MI_zero_idx[i]] = binomial_rng(N_counts_ppr[MI_zero_idx[i]], prob_seq_MI[MI_zero_idx[i]]);
    }else{
      N_counts_ppr[MI_zero_idx[i]] = trunc_binomial_rng(N_obs_max, prob_seq_1[MI_zero_idx[i]]);
      D_counts_ppr[MI_zero_idx[i]] = 0;
    }
  }
  for ( i in 1:N_ind_with_MI_nnzero){
    ind_log_prob_mi[MI_nnzero_idx[i]] = 0;
    N_counts_ppr[MI_nnzero_idx[i]] = trunc_binomial_rng(N_obs_max, prob_seq_any[MI_nnzero_idx[i]]);
    D_counts_ppr[MI_nnzero_idx[i]] = binomial_rng(N_counts_ppr[MI_nnzero_idx[i]], prob_seq_MI[MI_nnzero_idx[i]]);
  }
}

