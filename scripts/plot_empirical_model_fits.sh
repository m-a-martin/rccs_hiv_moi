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
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_age_sex_comm.Rds

#### ----------------------- ####
#### EMPIRICAL MI PREDICTORS ####
#### ----------------------- ####
Rscript scripts/plot_empirical_mi_predictors.R \
	--commTypeFit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_age_sex_comm.Rds \
	--sexpeverFit fit/211220_allreads_phsc_all_subgraphs_format_par_m_deep-phyloMI_sexpever_men.Rds \
	--varSelectFit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_var_select.Rds \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv  \
	--metadata output/211220_allreads_phsc_metadata.tsv \
	--out empirical_mi_predictors

# -------------------------------------------------------------------------------- #
#### +++++++++++++++++++++ ####
#### SUPPLEMENTARY FIGURES ####
#### +++++++++++++++++++++ ####
#### ---------------------- ####
#### EMPIRICAL AGE SEX COMM ####
#### ---------------------- ####
# TRACE
grep -v "prev. MIs" config/empirical_age_sex_comm_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_trace.R \
	--out empirical_age_sex_comm_trace \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_age_sex_comm.Rds \
	--plotParams config/tmp.tsv \
	--nRow 3
rm -rf config/tmp.tsv

# PAIRS
grep -v "prev. MIs" config/empirical_age_sex_comm_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_pairs.R \
	--out empirical_age_sex_comm_pairs \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_age_sex_comm.Rds \
	--plotParams config/tmp.tsv
rm -rf config/tmp.tsv

# PRIOR
grep -v "prev. MIs" config/empirical_age_sex_comm_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_fit_param_prior_posterior.R \
	--out empirical_age_sex_comm_prior \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_age_sex_comm.Rds  \
	--plotParams config/tmp.tsv \
	--nRowPer 3
rm -rf config/tmp.tsv

#### --------------------------------------------------------- ####
#### EMPIRICAL AGE, SEX, COMM_TYPE NO FREQ. CONTAMINANT FILTER ####
#### --------------------------------------------------------- ####
# TRACE
grep -v "prev. MIs" config/empirical_age_sex_comm_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_trace.R \
	--out empirical_age_sex_comm_incl_minor_trace \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_age_sex_comm_incl_minor.Rds \
	--plotParams config/tmp.tsv \
	--nRow 3
rm -rf config/tmp.tsv

# PAIRS
grep -v "prev. MIs" config/empirical_age_sex_comm_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_pairs.R \
	--out empirical_age_sex_comm_incl_minor_pairs \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_age_sex_comm_incl_minor.Rds \
	--plotParams config/tmp.tsv
rm -rf config/tmp.tsv

# PRIOR
grep -v "prev. MIs" config/empirical_age_sex_comm_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_fit_param_prior_posterior.R \
	--out empirical_age_sex_comm_incl_minor_prior \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_age_sex_comm_incl_minor.Rds  \
	--plotParams config/tmp.tsv \
	--nRowPer 3
rm -rf config/tmp.tsv

#### -------------------------------------------- ####
#### EMPIRICAL AGE SEX COMM  W/ SYN. AMPLICON DAT ####
#### -------------------------------------------- ####
# TRACE
grep -v "prev. MIs" config/empirical_age_sex_comm_plot_params.tsv | \
	grep -v "logit_prob_seq_ind_sd" \
	> config/tmp.tsv
Rscript scripts/plot_trace.R \
	--out empirical_age_sex_comm_amplicon_trace \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_amplicon_tmp_age_sex_comm_amplicon.Rds \
	--plotParams config/tmp.tsv \
	--nRow 3
rm -rf config/tmp.tsv

# PAIRS
grep -v "prev. MIs" config/empirical_age_sex_comm_plot_params.tsv | \
	grep -v "logit_prob_seq_ind_sd" \
	> config/tmp.tsv
Rscript scripts/plot_pairs.R \
	--out empirical_age_sex_comm_amplicon_pairs \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_amplicon_tmp_age_sex_comm_amplicon.Rds \
	--plotParams config/tmp.tsv
rm -rf config/tmp.tsv

# PRIOR
grep -v "prev. MIs" config/empirical_age_sex_comm_plot_params.tsv | \
	grep -v "logit_prob_seq_ind_sd" \
	> config/tmp.tsv
Rscript scripts/plot_fit_param_prior_posterior.R \
	--out empirical_age_sex_comm_amplicon_prior \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_amplicon_tmp_age_sex_comm_amplicon.Rds  \
	--plotParams config/tmp.tsv \
	--nRowPer 3
rm -rf config/tmp.tsv

#### ------------------------------------- ####
#### EMPIRICAL SEQUENCING TECHNOLOGY MODEL ####
#### ------------------------------------- ####
# TRACE
grep -v "prev. MIs" config/empirical_seq_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_trace.R \
	--out empirical_seqtech_trace \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_seq.Rds \
	--plotParams config/tmp.tsv \
	--nRow 3
rm -rf config/tmp.tsv

# PAIRS
grep -v "prev. MIs" config/empirical_seq_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_pairs.R \
	--out empirical_seqtech_pairs \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_seq.Rds \
	--plotParams config/tmp.tsv
rm -rf config/tmp.tsv

# PRIOR
grep -v "prev. MIs" config/empirical_seq_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_fit_param_prior_posterior.R \
	--out empirical_seqtech_prior \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_seq.Rds  \
	--plotParams config/tmp.tsv \
	--nRowPer 3
rm -rf config/tmp.tsv


#### --------------------------------------------------- ####
#### EMPIRICAL COMM_TYPE AND SEQUENCING_TECHNOLOGY MODEL ####
#### --------------------------------------------------- ####
# TRACE
grep -v "prev. MIs" config/empirical_comm_seq_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_trace.R \
	--out empirical_seqtech_comm_trace \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_comm_seq.Rds \
	--plotParams config/tmp.tsv \
	--nRow 3
rm -rf config/tmp.tsv

# PAIRS
grep -v "prev. MIs" config/empirical_comm_seq_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_pairs.R \
	--out empirical_seqtech_comm_pairs \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_comm_seq.Rds \
	--plotParams config/tmp.tsv
rm -rf config/tmp.tsv

# PRIOR
grep -v "prev. MIs" config/empirical_comm_seq_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_fit_param_prior_posterior.R \
	--out empirical_seqtech_comm_prior \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_comm_seq.Rds  \
	--plotParams config/tmp.tsv \
	--nRowPer 3
rm -rf config/tmp.tsv


#### ------------------------ ####
#### EMPIRICAL SEXPEVER MODEL ####
#### ------------------------ ####
# RESULTS ACROSS AGE GROUPS
Rscript scripts/plot_empirical_sexpever_pred.R \
	--sexpeverFit fit/211220_allreads_phsc_all_subgraphs_format_par_m_deep-phyloMI_sexpever_men.Rds \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv  \
	--out empirical_sexpever_men_prob_mi

# TRACE
grep -v "prev. MIs" config/empirical_sexpever_men_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_trace.R \
	--out empirical_sexpever_men_trace \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_m_deep-phyloMI_sexpever_men.Rds \
	--plotParams config/tmp.tsv \
	--nRow 4
rm -rf config/tmp.tsv

# PAIRS
grep -v "prev. MIs" config/empirical_sexpever_men_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_pairs.R \
	--out empirical_sexpever_men_pairs \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_m_deep-phyloMI_sexpever_men.Rds \
	--plotParams config/tmp.tsv
rm -rf config/tmp.tsv

# PRIOR
grep -v "prev. MIs" config/empirical_sexpever_men_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_fit_param_prior_posterior.R \
	--out empirical_sexpever_men_prior \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_m_deep-phyloMI_sexpever_men.Rds   \
	--plotParams config/tmp.tsv \
	--nRowPer 3
rm -rf config/tmp.tsv


#### ----------------------------------------- ####
#### 4. EXTENDED MODEL WITH VARIABLE SELECTION ####
#### ----------------------------------------- ####
# TRACE
grep -v "prev. MIs" config/empirical_var_select_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_trace.R \
	--out empirical_var_select_trace \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_var_select.Rds \
	--plotParams config/tmp.tsv \
	--nRow 10
rm -rf config/tmp.tsv


# PAIRS
grep -v "prev. MIs" config/empirical_var_select_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_pairs.R \
	--out empirical_var_select_pairs \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_var_select.Rds \
	--plotParams config/tmp.tsv
rm -rf config/tmp.tsv


# PRIOR 
grep -v "prev. MIs" config/empirical_var_select_plot_params.tsv \
	> config/tmp.tsv
Rscript scripts/plot_fit_param_prior_posterior.R \
	--out empirical_var_select_prior \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_var_select.Rds  \
	--plotParams config/tmp.tsv \
	--nRowPer 10
rm -rf config/tmp.tsv

