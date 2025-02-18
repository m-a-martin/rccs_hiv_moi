#!/usr/bin/env bash

#### ------------------------------------------------------------ ####
#### 1. EXTENDED MODEL WITH PROB. MI BASED ON AGE, SEX, COMM_TYPE ####
#### ------------------------------------------------------------ ####
Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--inputType full \
	--stan stan/deep-phyloMI.stan \
	--filter "id_subgraph_reads > 0 & window_type == 'unique'" \
	--miDesignMatrix age_cat_coarse sex comm_type \
	--seqDesignMatrix "sequencing_technology" "log10_copies_shared" "sequencing_technology*log10_copies" \
	--scaleVars "log10_copies_shared" "log10_copies:sequencing_technology" \
	--strataPrevs "comm_type == 'fishing'" "comm_type == 'inland'" \
	--strataPrevRatios "comm_type == 'fishing':comm_type=='inland'" \
	--outAppend age_sex_comm


# finally, a bit of post-hoc analyses
# poststratification weighting
Rscript scripts/poststratification.R \
	--metadata output/211220_allreads_phsc_metadata.tsv \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_age_sex_comm.Rds \
	--strataVars age_cat_coarse sex comm_type \
	--strataPrevs "comm_type == 'fishing'" "comm_type == 'inland'" \
	--strataPrevRatios "comm_type == 'fishing':comm_type=='inland'"

# MI classification threshold
Rscript scripts/calc_n_mi_threshold.R \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_age_sex_comm.Rds \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv


#### ---------------------------------------------------------------------------------------- ####
#### 2. EXTENDED MODEL WITH PROB. MI BASED ON AGE, SEX, COMM_TYPE NO FREQ. CONTAMINANT FILTER ####
#### ---------------------------------------------------------------------------------------- ####
Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--inputType full \
	--filter "window_type == 'unique' & id_subgraph_reads_all > 0 & subgraph_reads_all >= 3" \
	--stan stan/deep-phyloMI.stan \
	--miDesignMatrix age_cat_coarse sex comm_type \
	--seqDesignMatrix "sequencing_technology" "log10_copies_shared" "sequencing_technology*log10_copies" \
	--scaleVars "log10_copies_shared" "log10_copies:sequencing_technology" \
	--outAppend age_sex_comm_incl_minor



#### --------------------------------------------------------------------------------- ####
#### 3. EXTENDED MODEL WITH PROB. MI BASED ON AGE, SEX, COMM_TYPE W/ SYN. AMPLICON DAT ####
#### --------------------------------------------------------------------------------- ####
sed 's/target += cauchy_lpdf(logit_prob_seq_ind_sd/\/\/ target += cauchy_lpdf(logit_prob_seq_ind_sd/g' \
	stan/deep-phyloMI.stan |
	sed 's/real<lower = 0> logit_prob_seq_ind_sd/\/\/ real<lower = 0> logit_prob_seq_ind_sd/g' | \
	sed 's/0, logit_prob_seq_ind_sd/0, 0.7/g' |
	> stan/tmp.stan

Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par_amplicon.tsv \
	--inputType window \
	--filter "TRUE" \
	--stan stan/tmp.stan \
	--miDesignMatrix age_cat_coarse sex comm_type \
	--seqDesignMatrix "sequencing_technology" "log10_copies_shared"  "sequencing_technology*log10_copies" \
	--scaleVars "log10_copies_shared" "log10_copies:sequencing_technology" \
	--outAppend age_sex_comm_amplicon

rm -rf stan/tmp.stan
#### -------------------------------------------------------------- ####
#### 4. EXTENDED MODEL WITH PROB. MI BASED ON SEQUENCING TECHNOLOGY ####
#### ------------------------------------------------------------- ####
Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--inputType full \
	--stan stan/deep-phyloMI.stan \
	--filter "id_subgraph_reads > 0 & window_type == 'unique'" \
	--miDesignMatrix sequencing_technology \
	--seqDesignMatrix "sequencing_technology" "log10_copies_shared"  "sequencing_technology*log10_copies" \
	--scaleVars '"log10_copies_shared" log10_copies:sequencing_technology' \
	--multivariateRiskRatios "sequencing_technology=='bait_capture':sequencing_technology=='amplicon'" \
	--outAppend seq


#### ---------------------------------------------------------------------------- ####
#### 5. EXTENDED MODEL WITH PROB. MI BASED ON COMM_TYPE AND SEQUENCING_TECHNOLOGY ####
#### ---------------------------------------------------------------------------- ####
Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--inputType full \
	--stan stan/deep-phyloMI.stan \
	--filter "id_subgraph_reads > 0 & window_type == 'unique'" \
	--miDesignMatrix comm_type sequencing_technology \
	--seqDesignMatrix "sequencing_technology" "log10_copies_shared"  "sequencing_technology*log10_copies" \
	--scaleVars "log10_copies_shared" 'log10_copies:sequencing_technology' \
	--multivariateRiskRatios "sequencing_technology=='bait_capture':sequencing_technology=='amplicon'" \
	--outAppend comm_seq

#### ----------------------------------------------------------------- ####
#### 6. EXTENDED MODEL WITH PROB. MI DEPENDENT ON COMM_TYPE & SEXPEVER ####
#### ----------------------------------------------------------------- ####
# ONLY PERFORM THIS ANALYSIS AMONG MEN
Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par_m.tsv \
	--inputType full \
	--stan stan/deep-phyloMI.stan \
	--filter "id_subgraph_reads > 0 & window_type == 'unique'" \
	--seqDesignMatrix "sequencing_technology" "log10_copies_shared"  "sequencing_technology*log10_copies" \
	--miDesignMatrix "comm_type" "plhiv_sexpever_std" "comm_type*plhiv_sexpever_std" \
	--scaleVars "log10_copies_shared" 'log10_copies:sequencing_technology' \
	--miMissingDat \
		'plhiv_sexpever_std:93;3;60;plhiv_sexpeverMany_shape;plhiv_sexpeverMany_scale;plhiv_sexpever_mean' \
		'plhiv_sexpever_std:100;1;60;plhiv_sexpeverMissing_shape;plhiv_sexpeverMissing_scale;plhiv_sexpever_mean' \
	--outAppend sexpever_men

Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par_m_complete.tsv \
	--inputType full \
	--stan stan/deep-phyloMI.stan \
	--filter "id_subgraph_reads > 0 & window_type == 'unique'" \
	--seqDesignMatrix "sequencing_technology" "log10_copies_shared"  "sequencing_technology*log10_copies" \
	--miDesignMatrix "comm_type" "plhiv_sexpever_std" "comm_type*plhiv_sexpever_std" \
	--scaleVars "log10_copies_shared" 'log10_copies:sequencing_technology' \
	--outAppend sexpever_men_complete


#### ----------------------------------------- ####
#### 7. EXTENDED MODEL WITH VARIABLE SELECTION ####
#### ----------------------------------------- ####
Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--inputType full \
	--stan stan/deep-phyloMI.stan \
	--miCoeffPriors shrinkage \
	--filter "id_subgraph_reads > 0 & window_type == 'unique'" \
	--miDesignMatrix \
		"sequencing_technology" \
		"round" \
		"male_circumcision" \
		"comm_type" \
		"sex" \
		"age_cat_coarse" \
		"in_migrant" \
		"barworker" \
		"sex*comm_type"  \
		"age_cat_coarse*comm_type" \
		"in_migrant*comm_type" \
		"barworker*comm_type" \
	--seqDesignMatrix "sequencing_technology" "log10_copies_shared"  "sequencing_technology*log10_copies" \
	--scaleVars "log10_copies_shared" "log10_copies:sequencing_technology" \
	--strataPrevRatios "comm_type=='inland'&in_migrant==TRUE:comm_type=='inland'&in_migrant==FALSE"  "comm_type=='fishing'&sex=='M':comm_type=='fishing'&sex=='F'" \
	--multivariateRiskRatios "comm_type=='fishing':comm_type=='inland'" "comm_type=='inland'&in_migrant==TRUE:comm_type=='inland'&in_migrant==FALSE"\
	--outAppend var_select


#### -------------------------------- ####
#### 8. FORMAT PARAMETER CONFIG FILES ####
#### -------------------------------- ####
# AGE SEX COMM
bash scripts/format_params.sh \
	config/empirical_extended_plot_params_base.tsv \
	age_sex_comm

# SEQ
bash scripts/format_params.sh \
	config/empirical_extended_plot_params_base.tsv \
	seq

# COMM SEQ
bash scripts/format_params.sh \
	config/empirical_extended_plot_params_base.tsv \
	comm_seq

# SEXPEVER MEN
bash scripts/format_params.sh \
	config/empirical_extended_plot_params_base.tsv \
	sexpever_men \
	m_

# VAR SELECT
bash scripts/format_params.sh \
	config/empirical_extended_plot_params_base.tsv \
	var_select

