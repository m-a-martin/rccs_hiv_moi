#!/usr/bin/env bash

#### ++++++++++++++++++ ####
#### MAIN PAPER FIGURES ####
#### ++++++++++++++++++ ####
#### ------------------------ ####
#### EMPIRICAL FULL MODEL FIT ####
#### ------------------------ ####
Rscript scripts/plot_empirical_fit.R \
	--out empirical_full_fit \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_full_model_.Rds

#### ----------------------- ####
#### EMPIRICAL MI PREDICTORS ####
#### ----------------------- ####
Rscript scripts/plot_empirical_mi_predictors.R \
	--commTypeFit fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_comm_type.Rds \
	--sexpeverFit fit/211220_allreads_phsc_all_subgraphs_format_par_m_extended_model_sexpever_men.Rds \
	--varSelectFit fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_hsp_var_select.Rds \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv  \
	--out empirical_mi_predictors

# -------------------------------------------------------------------------------- #
#### +++++++++++++++++++++ ####
#### SUPPLEMENTARY FIGURES ####
#### +++++++++++++++++++++ ####
#### -------------------- ####
#### EMPIRICAL FULL MODEL ####
#### -------------------- ####
# TRACE
grep -v "prev. MIs" config/full_empirical_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_trace.R \
	--out empirical_full_trace \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_full_model_.Rds \
	--plotParams config/tmp.tsv \
	--nRow 3
rm -rf config/tmp.tsv

# PAIRS
grep -v "prev. MIs" config/full_empirical_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_pairs.R \
	--out empirical_full_pairs \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_full_model_.Rds \
	--plotParams config/tmp.tsv
rm -rf config/tmp.tsv

# PRIOR
grep -v "prev. MIs" config/full_empirical_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_fit_param_prior_posterior.R \
	--out empirical_full_prior \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_full_model_.Rds  \
	--plotParams config/tmp.tsv \
	--nRowPer 3
rm -rf config/tmp.tsv

# PARAMS
grep -v "prev. MIs" config/full_empirical_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_fit_params.R \
	--out empirical_full_params \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_full_model_.Rds \
	--plotParams config/tmp.tsv \
	--nRowPer 3
rm -rf config/tmp.tsv


#### ------------------------- ####
#### EMPIRICAL SEQ_TECH MODEL ####
#### ------------------------- ####
# TRACE
grep -v "prev. MIs" config/seqtech_empirical_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_trace.R \
	--out empirical_seqtech_trace \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_sequencing_technology.Rds \
	--plotParams config/tmp.tsv \
	--nRow 3
rm -rf config/tmp.tsv

# PAIRS
grep -v "prev. MIs" config/seqtech_empirical_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_pairs.R \
	--out empirical_seqtech_pairs \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_sequencing_technology.Rds \
	--plotParams config/tmp.tsv
rm -rf config/tmp.tsv

# PRIOR
grep -v "prev. MIs" config/seqtech_empirical_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_fit_param_prior_posterior.R \
	--out empirical_seqtech_prior \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_sequencing_technology.Rds  \
	--plotParams config/tmp.tsv \
	--nRowPer 3
rm -rf config/tmp.tsv

# PARAMS
grep -v "prev. MIs" config/seqtech_empirical_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_fit_params.R \
	--out empirical_seqtech_params \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_sequencing_technology.Rds \
	--plotParams config/tmp.tsv \
	--nRowPer 3
rm -rf config/tmp.tsv


#### ---------------------------------- ####
#### EMPIRICAL SEQ_TECH_COMM_TYPE MODEL ####
#### ---------------------------------- ####
# TRACE
grep -v "prev. MIs" config/sequencing_technology_comm_type_empirical_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_trace.R \
	--out empirical_seqtech_comm_trace \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_sequencing_technology_comm_type.Rds \
	--plotParams config/tmp.tsv \
	--nRow 3
rm -rf config/tmp.tsv

# PAIRS
grep -v "prev. MIs" config/sequencing_technology_comm_type_empirical_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_pairs.R \
	--out empirical_seqtech_comm_pairs \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_sequencing_technology_comm_type.Rds \
	--plotParams config/tmp.tsv
rm -rf config/tmp.tsv

# PRIOR
grep -v "prev. MIs" config/sequencing_technology_comm_type_empirical_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_fit_param_prior_posterior.R \
	--out empirical_seqtech_comm_prior \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_sequencing_technology_comm_type.Rds  \
	--plotParams config/tmp.tsv \
	--nRowPer 3
rm -rf config/tmp.tsv

# PARAMS
grep -v "prev. MIs" config/sequencing_technology_comm_type_empirical_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_fit_params.R \
	--out empirical_seqtech_comm_params \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_sequencing_technology_comm_type.Rds \
	--plotParams config/tmp.tsv \
	--nRowPer 3
rm -rf config/tmp.tsv


#### ------------------------- ####
#### EMPIRICAL COMM_TYPE MODEL ####
#### ------------------------- ####
# TRACE
grep -v "prev. MIs" config/commtype_empirical_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_trace.R \
	--out empirical_commtype_trace \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_comm_type.Rds \
	--plotParams config/tmp.tsv \
	--nRow 3
rm -rf config/tmp.tsv

# PAIRS
grep -v "prev. MIs" config/commtype_empirical_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_pairs.R \
	--out empirical_commtype_pairs \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_comm_type.Rds \
	--plotParams config/tmp.tsv
rm -rf config/tmp.tsv

# PRIORS
grep -v "prev. MIs" config/commtype_empirical_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_fit_param_prior_posterior.R \
	--out empirical_commtype_prior \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_comm_type.Rds  \
	--plotParams config/tmp.tsv \
	--nRowPer 3
rm -rf config/tmp.tsv

# PARAMS
grep -v "prev. MIs" config/commtype_empirical_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_fit_params.R \
	--out empirical_commtype_params \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_comm_type.Rds \
	--plotParams config/tmp.tsv \
	--nRowPer 3
rm -rf config/tmp.tsv


#### ------------------------ ####
#### EMPIRICAL SEXPEVER MODEL ####
#### ------------------------ ####
# IMPUTED VALUES
Rscript scripts/plot_empirical_imputed_sexpever.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par_m.tsv \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_m_extended_model_sexpever_men.Rds \
	--out empirical_sexpever_men_imputed

# TRACE
Rscript scripts/plot_trace.R \
	--out empirical_sexpever_men_trace \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_m_extended_model_sexpever_men.Rds \
	--plotParams config/sexpever_men_empirical_plot_params.tsv \
	--nRow 3

# PAIRS
Rscript scripts/plot_pairs.R \
	--out empirical_sexpever_men_pairs \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_m_extended_model_sexpever_men.Rds \
	--plotParams config/sexpever_men_empirical_plot_params.tsv

# PRIOR
Rscript scripts/plot_fit_param_prior_posterior.R \
	--out empirical_sexpever_men_prior \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_m_extended_model_sexpever_men.Rds  \
	--plotParams config/sexpever_men_empirical_plot_params.tsv \
	--nRowPer 3

# PARAMS
Rscript scripts/plot_fit_params.R \
	--out empirical_sexpever_men_params \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_m_extended_model_sexpever_men.Rds \
	--plotParams config/sexpever_men_empirical_plot_params.tsv \
	--nRowPer 3

# PREDICTED PROB. MI ACROSS AGE CATEGORIES
Rscript scripts/plot_empirical_sexpever_pred.R \
	--sexpeverFit fit/211220_allreads_phsc_all_subgraphs_format_par_m_extended_model_sexpever_men.Rds \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv  \
	--out empirical_sexpever_men_probMI


#### ----------------------------------------- ####
#### 4. EXTENDED MODEL WITH VARIABLE SELECTION ####
#### ----------------------------------------- ####
# TRACE
Rscript scripts/plot_trace.R \
	--out empirical_var_select_trace \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_hsp_var_select.Rds \
	--plotParams config/var_select_empirical_plot_params.tsv \
	--nRow 6

# PAIRS
Rscript scripts/plot_pairs.R \
	--out empirical_var_select_pairs \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_hsp_var_select.Rds \
	--plotParams config/var_select_empirical_plot_params.tsv

# PRIOR 
Rscript scripts/plot_fit_param_prior_posterior.R \
	--out empirical_var_select_prior \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_hsp_var_select.Rds  \
	--plotParams config/var_select_empirical_plot_params.tsv \
	--nRowPer 6

# PARAMS
# some params plotted in main figure, so excluding
grep -v "w\\[" config/var_select_empirical_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_fit_params.R \
	--out empirical_var_select_params \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_hsp_var_select.Rds  \
	--plotParams config/tmp.tsv \
	--nRowPer 3
rm -rf config/tmp.tsv
