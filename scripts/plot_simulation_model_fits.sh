#!/usr/bin/env bash

#### ++++++++++++++++++ ####
#### MAIN PAPER FIGURES ####
#### ++++++++++++++++++ ####
#### ------------------- ####
#### FULL SIMULATION FIT ####
#### ------------------- ####
grep -v "delta_0" config/full_simulation_plot_params.tsv |\
	grep -v "logit_prob_MI_fnr" | \
	grep -v "logit_prob_MI_fpr" | \
	sed "s/VL/viral load/g" \
	> config/tmp.tsv
Rscript scripts/plot_simulation_fit.R \
	--out full_simulation_full_fit \
	--dat simulations/full_simulation.tsv \
	--params simulations/full_simulation_params.tsv \
	--fit fit/full_simulation_full_model_.Rds \
	--plotParams config/tmp.tsv
rm -rf config/tmp.tsv


# -------------------------------------------------------------------------------- #
#### +++++++++++++++++++++ ####
#### SUPPLEMENTARY FIGURES ####
#### +++++++++++++++++++++ ####
#### -------------------------- ####
#### BASE SIMULATION BASE MODEL ####
#### -------------------------- ####
# TRACE
grep -v "prev. MIs" config/base_simulation_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_trace.R \
	--out base_simulation_base_trace \
	--fit fit/base_simulation_base_model_.Rds \
	--plotParams config/tmp.tsv \
	--nRow 3
rm -rf config/tmp.tsv

# PAIRS
grep -v "prev. MIs" config/base_simulation_plot_params.tsv > config/tmp.tsv
Rscript scripts/plot_pairs.R \
	--out base_simulation_base_pairs \
	--fit fit/base_simulation_base_model_.Rds \
	--plotParams config/tmp.tsv 
rm -rf config/tmp.tsv

# PRIOR
grep -v "prev. MIs" config/base_simulation_plot_params.tsv > config/tmp.tsv
Rscript scripts/plot_fit_param_prior_posterior.R \
	--out base_simulation_base_prior \
	--fit fit/base_simulation_base_model_.Rds \
	--plotParams config/tmp.tsv	
rm -rf config/tmp.tsv

# FIT
grep -v "delta_0" config/base_simulation_plot_params.tsv | \
	sed 's/VL/viral load/g' \
	> config/tmp.tsv
Rscript scripts/plot_simulation_fit.R \
	--out base_simulation_base_fit \
	--dat simulations/base_simulation.tsv \
	--params simulations/base_simulation_params.tsv \
	--fit fit/base_simulation_base_model_.Rds \
	--plotParams config/tmp.tsv \
	--ncolParams 4
rm -rf config/tmp.tsv


#### -------------------------- ####
#### FULL_SIMULATION BASE MODEL ####
#### -------------------------- ####
# TRACE
grep -v "prev. MIs" config/base_simulation_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_trace.R \
	--out full_simulation_base_trace \
	--fit fit/full_simulation_base_model_.Rds \
	--plotParams config/tmp.tsv \
	--nRow 3
rm -rf config/tmp.tsv

# PAIRS
grep -v "prev. MIs" config/base_simulation_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_pairs.R \
	--out full_simulation_base_pairs \
	--fit fit/full_simulation_base_model_.Rds \
	--plotParams config/tmp.tsv 
rm -rf config/tmp.tsv

# PRIOR
grep -v "prev. MIs" config/base_simulation_plot_params.tsv > config/tmp.tsv
Rscript scripts/plot_fit_param_prior_posterior.R \
	--out full_simulation_base_prior \
	--fit fit/full_simulation_base_model_.Rds \
	--plotParams config/tmp.tsv	
rm -rf config/tmp.tsv

# PARAMS 
grep -v "delta_0" config/base_simulation_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_fit_params.R \
	--out full_simulation_base_params \
	--params simulations/full_simulation_params.tsv \
	--fit fit/full_simulation_base_model_.Rds \
	--plotParams config/tmp.tsv \
	--nRowPer 1
rm -rf config/tmp.tsv


#### -------------------------- ####
#### FULL_SIMULATION FULL MODEL ####
#### -------------------------- ####
# TRACE
grep -v "prev. MIs" config/full_simulation_plot_params.tsv | \
	grep -v "logit_prob_MI_fnr" | \
	grep -v "logit_prob_MI_fpr" \
	> config/tmp.tsv
Rscript scripts/plot_trace.R \
	--out full_simulation_full_trace \
	--fit fit/full_simulation_full_model_.Rds \
	--plotParams config/tmp.tsv \
	--nRow 3
rm -rf config/tmp.tsv

# PAIRS
grep -v "prev. MIs" config/full_simulation_plot_params.tsv  | \
	grep -v "logit_prob_MI_fnr" | \
	grep -v "logit_prob_MI_fpr" \
	> config/tmp.tsv
Rscript scripts/plot_pairs.R \
	--out full_simulation_full_pairs \
	--fit fit/full_simulation_full_model_.Rds \
	--plotParams config/tmp.tsv 
rm -rf config/tmp.tsv

# PRIOR
grep -v "prev. MIs" config/full_simulation_plot_params.tsv  | \
	grep -v "logit_prob_MI_fnr" | \
	grep -v "logit_prob_MI_fpr" \
	> config/tmp.tsv
Rscript scripts/plot_fit_param_prior_posterior.R \
	--out full_simulation_full_prior \
	--fit fit/full_simulation_full_model_.Rds \
	--plotParams config/tmp.tsv	\
	--nRowPer 2
rm -rf config/tmp.tsv

# PARAMS 
grep -v "prev. MIs" config/full_simulation_plot_params.tsv  | \
	grep -v "logit_prob_MI_fnr" | \
	grep -v "logit_prob_MI_fpr" \
	> config/tmp.tsv
Rscript scripts/plot_fit_params.R \
	--out full_simulation_full_params \
	--params simulations/full_simulation_params.tsv \
	--fit fit/full_simulation_full_model_.Rds \
	--plotParams config/tmp.tsv \
	--nRowPer 1
rm -rf config/tmp.tsv

#### --------------------------------------------- ####
####  FULL_SIMULATION FULL MODEL DELTA SENSITIVITY ####
#### --------------------------------------------- ####
grep -v "logit_prob_MI" config/full_simulation_sensitivity_plot_params.tsv |\
	grep -v "logit_prob_MI_fnr" | \
	grep -v "logit_prob_MI_fpr" \
	> config/tmp.tsv
Rscript scripts/plot_fit_params.R \
	--out full_simulation_sensitivity_delta \
	--fit 'fit/sensitivity/*delta*.Rds' \
	--params 'simulations/sensitivity/*delta*_params.tsv' \
	--plotParams config/tmp.tsv
rm -rf config/tmp.tsv


#### --------------------------------------------- ####
#### FULL_SIMULATION FULL MODEL LAMBDA SENSITIVITY ####
#### --------------------------------------------- ####
grep -v "logit_prob_MI" config/full_simulation_sensitivity_plot_params.tsv |\
	grep -v "logit_prob_MI_fnr" | \
	grep -v "logit_prob_MI_fpr" \
	> config/tmp.tsv
Rscript scripts/plot_fit_params.R \
	--out full_simulation_sensitivity_lambda \
	--fit 'fit/sensitivity/*fnr*.Rds' \
	--params 'simulations/sensitivity/*fnr*_params.tsv' \
	--plotParams config/tmp.tsv
rm -rf config/tmp.tsv


#### ---------------------------------------------- ####
#### FULL_SIMULATION FULL MODEL EPSILON SENSITIVITY ####
#### ---------------------------------------------- ####
grep -v "logit_prob_MI" config/full_simulation_sensitivity_plot_params.tsv |\
	grep -v "logit_prob_MI_fnr" | \
	grep -v "logit_prob_MI_fpr" \
	> config/tmp.tsv
Rscript scripts/plot_fit_params.R \
	--out full_simulation_sensitivity_epsilon \
	--fit 'fit/sensitivity/*fpr*.Rds' \
	--params 'simulations/sensitivity/*fpr*_params.tsv' \
	--plotParams config/tmp.tsv
rm -rf config/tmp.tsv


#### ---------------------------------------- ####
#### EXTENDED SIMULATION EXTENDED MODEL TRACE ####
#### ---------------------------------------- ####
# TRACE
grep -v "prev. MIs" config/extended_simulation_plot_params.tsv | \
	grep -v "logit_prob_MI_fnr" | \
	grep -v "logit_prob_MI_fpr" \
	> config/tmp.tsv
Rscript scripts/plot_trace.R \
	--out extended_simulation_extended_trace \
	--fit fit/extended_simulation_extended_model_hsp_.Rds \
	--plotParams config/tmp.tsv \
	--nRow 4
rm -rf config/tmp.tsv


# PAIRS
grep -v "prev. MIs" config/extended_simulation_plot_params.tsv | \
	grep -v "logit_prob_MI_fnr" | \
	grep -v "logit_prob_MI_fpr" \
	> config/tmp.tsv
Rscript scripts/plot_pairs.R \
	--out extended_simulation_extended_pairs \
	--fit fit/extended_simulation_extended_model_hsp_.Rds \
	--plotParams config/tmp.tsv
rm -rf config/tmp.tsv

# PRIOR
grep -v "prev. MIs" config/extended_simulation_plot_params.tsv  | \
	grep -v "logit_prob_MI_fnr" | \
	grep -v "logit_prob_MI_fpr" \
	> config/tmp.tsv
Rscript scripts/plot_fit_param_prior_posterior.R \
	--out extended_simulation_extended_prior \
	--fit fit/extended_simulation_extended_model_hsp_.Rds \
	--plotParams config/tmp.tsv	\
	--nRowPer 2
rm -rf config/tmp.tsv

# FIT
grep -v "prev. MIs" config/extended_simulation_plot_params.tsv | \
	grep -v "logit_prob_MI_fnr" | \
	grep -v "logit_prob_MI_fpr" \
	> config/tmp.tsv
Rscript scripts/plot_simulation_fit.R \
	--out extended_simulation_extended_fit \
	--fit fit/extended_simulation_extended_model_hsp_.Rds \
	--dat simulations/extended_simulation.tsv \
	--params simulations/extended_simulation_params.tsv \
	--plotParams config/tmp.tsv \
	--ncolParams 5
rm -rf config/tmp.tsv

