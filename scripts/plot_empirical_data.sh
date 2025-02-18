#!/usr/bin/env bash

#### ++++++++++++++++++ ####
#### MAIN PAPER FIGURES ####
#### ++++++++++++++++++ ####

#### ---------------------- ####
#### EMPIRICAL DATA SUMMARY ####
#### ---------------------- ####
# empirical data summary
Rscript scripts/plot_empirical_data_sum.R \
	--out empirical_data_sum \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--filter "id_subgraph_reads > 0 & window_type == 'unique'" \
	--repTree1 output/trees/AID2211-fq1_1800_to_2049.Rda \
	--repTree2 output/trees/AID2642-fq2_1800_to_2049.Rda


# -------------------------------------------------------------------------------- #
#### +++++++++++++++++++++ ####
#### SUPPLEMENTARY FIGURES ####
#### +++++++++++++++++++++ ####
#### ----------------- ####
#### SEXPEVER OBSERVED ####
#### ----------------- ####
Rscript scripts/plot_mean_obs_sexpever.R \
	--dat 'data/input_metadata.tsv' \
	--out 'mean_obs_sexpever'

#### -------------- ####
#### VL SEQ SUCCESS ####
#### -------------- ####
# sequencing success as a function of viral load
Rscript scripts/plot_vl_seq_success.R \
	--out vl_seq_success \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--filter "id_subgraph_reads > 0 & window_type == 'unique'"


#### ------------------- ####
#### SEXPEVER STD. CURVE ####
#### ------------------- ####
# includes naive imputation of missing values
Rscript scripts/plot_sexpever_std_curve.R \
	--dat output/211220_allreads_phsc_metadata.tsv \
	--out sexpever_std_curve




