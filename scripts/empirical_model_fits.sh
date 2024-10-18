#!/usr/bin/env bash


#### ------------------------------------- ####
#### 1. FULL MODEL WITH IDENTICAL PROB. MI ####
#### ------------------------------------- ####
Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--inputType full \
	--stan stan/full_model.stan \
	--seqDesignMatrix "sequencing_technology" "sequencing_technology:log10_copies" \
	--scaleVars "log10_copies:sequencing_technology" \
	--designRisks "logit_prob_MI"

# need to format plot params file based on design matrix
n_intercept=$(grep -v "\\:" fit/211220_allreads_phsc_all_subgraphs_format_par_full_model__design_cols.tsv | \
	grep -v "col" | wc -l)
cat \
	config/full_empirical_plot_params_base.tsv \
	<(echo "\n") \
	<(awk -F'\t' '{print "logit_prob_seq_coeffs["NR"]\texpression(atop(\x27logit(seq. success):\x27, paste(\x27"$1" coeff. (\x27,alpha^\"" $1 "\", \x27)\x27)))\talpha^"$1}' \
		<(awk -F'\t' '{if ($2=="seq") print $1}' \
			<(grep -v "\\:" fit/211220_allreads_phsc_all_subgraphs_format_par_full_model__design_cols.tsv))) \
	<(awk -F'\t' -v n=$n_intercept '{print "logit_prob_seq_coeffs["NR+n"]\texpression(atop(\x27logit(seq. success):\x27, paste(\x27"$1" coeff. (\x27,alpha^\"" $1 "\", \x27)\x27)))\talpha^"$1}' \
		<(awk -F'\t' '{if ($2=="seq") print $1}' \
			<(grep "\\:" fit/211220_allreads_phsc_all_subgraphs_format_par_full_model__design_cols.tsv))) |\
	sed '/^$/d' | \
	sed 's/sequencing_technology//g' |	\
	sed 's/bait_capture/bait capture/g' | \
	sed 's/:log10_copies/-vl/g' | \
	sed 's/-vl coeff/ VL coeff/g' | \
	sed 's/alpha^\"bait capture/alpha^\"bait/g' | \
	sed 's/alpha^bait capture/alpha^bait/g' \
	> config/tmp.tsv

# sort rows
cat \
	<(head -n 1 config/tmp.tsv) \
	<(tail -n +2 config/tmp.tsv | tail -n +2 config/tmp.tsv | grep "logit_prob_seq" | sort) \
	<(grep lambda config/tmp.tsv) \
	<(grep epsilon config/tmp.tsv) \
	<(grep prob_MI config/tmp.tsv | grep -v epsilon | grep -v lambda) \
	> config/full_empirical_plot_params.tsv
rm -rf config/tmp.tsv

# finally, a bit of post-hoc analyses
Rscript scripts/calc_n_mi_threshold.R \
	--fit fit/211220_allreads_phsc_all_subgraphs_format_par_full_model_.Rds \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv


#### ------------------------------------------------------------ ####
#### 2. EXTENDED MODEL WITH PROB. MI DEPENDENT ON SEQUENCING TECH ####
#### ------------------------------------------------------------ ####
Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--inputType full \
	--stan stan/extended_model.stan \
	--miDesignMatrix sequencing_technology \
	--seqDesignMatrix "sequencing_technology" "sequencing_technology:log10_copies" \
	--scaleVars 'log10_copies:sequencing_technology' \
	--designRisks "logit_prob_MI_coeffs[1]+logit_prob_MI" "logit_prob_MI" \
	--designRiskRatios "logit_prob_MI_coeffs[1]+logit_prob_MI:logit_prob_MI" \
	--outAppend sequencing_technology

# need to format plot params file based on design matrix
n_intercept=$(
	awk -F'\t' '{if ($2=="seq") print $0}' fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_sequencing_technology_design_cols.tsv | \
		grep -v "\\:"  | \
		grep -v "col" | wc -l)

cat \
	config/extended_empirical_plot_params_base.tsv \
	<(echo "\n") \
	<(awk -F'\t' '{print "logit_prob_seq_coeffs["NR"]\texpression(atop(\x27logit(seq. success):\x27, paste(\x27"$1" coeff. (\x27,alpha^\"" $1 "\", \x27)\x27)))\talpha^"$1}' \
		<(awk -F'\t' '{if ($2=="seq") print $1}' \
			<(grep -v "\\:" fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_sequencing_technology_design_cols.tsv))) \
	<(awk -F'\t' -v n=$n_intercept '{print "logit_prob_seq_coeffs["NR+n"]\texpression(atop(\x27logit(seq. success):\x27, paste(\x27"$1" coeff. (\x27,alpha^\"" $1 "\", \x27)\x27)))\talpha^"$1}' \
		<(awk -F'\t' '{if ($2=="seq") print $1}' \
			<(grep "\\:" fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_sequencing_technology_design_cols.tsv))) \
	<(awk -F'\t' '{print "logit_prob_MI_coeffs["NR"]\t\"expression(atop(\x27logit(prob. MI):\x27, paste(\x27"$1" coeff. (\x27,beta["NR"], \x27)\x27)))\"\tbeta_"NR}' \
		<(awk -F'\t' '{if ($2=="mi") print $1}' fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_sequencing_technology_design_cols.tsv)) |\
	sed '/^$/d' | \
	sed 's/sequencing_technology//g' |	\
	sed 's/bait_capture/bait capture/g' | \
	sed 's/:log10_copies/-vl/g' | \
	sed 's/-vl coeff/ VL coeff/g' | \
	sed 's/alpha^\"bait capture/alpha^\"bait/g' | \
	sed 's/alpha^bait capture/alpha^bait/g' \
	> config/tmp.tsv

# sort rows
cat \
	<(head -n 1 config/tmp.tsv) \
	<(tail -n +2 config/tmp.tsv | tail -n +2 config/tmp.tsv | grep "logit_prob_seq" | sort) \
	<(grep lambda config/tmp.tsv) \
	<(grep epsilon config/tmp.tsv) \
	<(grep prob_MI config/tmp.tsv | grep -v epsilon | grep -v lambda) \
	> config/seqtech_empirical_plot_params.tsv
rm -rf config/tmp.tsv


#### ------------------------------------------------------------------------ ####
#### 3. EXTENDED MODEL WITH PROB. MI DEPENDENT ON COMM_TYPE & SEQUENCING_TECH ####
#### ------------------------------------------------------------------------ ####
Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--inputType full \
	--stan stan/extended_model.stan \
	--miDesignMatrix comm_type sequencing_technology \
	--seqDesignMatrix "sequencing_technology" "sequencing_technology:log10_copies" \
	--scaleVars 'log10_copies:sequencing_technology' \
	--populationRiskRatios "sequencing_technology='bait_capture':sequencing_technology='amplicon'"
	--outAppend sequencing_technology_comm_type


# need to format plot params file based on design matrix
n_intercept=$(
	awk -F'\t' '{if ($2=="seq") print $0}' fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_sequencing_technology_comm_type_design_cols.tsv | \
		grep -v "\\:"  | \
		grep -v "col" | wc -l)

cat \
	config/extended_empirical_plot_params_base.tsv \
	<(echo "\n") \
	<(awk -F'\t' '{print "logit_prob_seq_coeffs["NR"]\texpression(atop(\x27logit(seq. success):\x27, paste(\x27"$1" coeff. (\x27,alpha^\"" $1 "\", \x27)\x27)))\talpha^"$1}' \
		<(awk -F'\t' '{if ($2=="seq") print $1}' \
			<(grep -v "\\:" fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_sequencing_technology_comm_type_design_cols.tsv))) \
	<(awk -F'\t' -v n=$n_intercept '{print "logit_prob_seq_coeffs["NR+n"]\texpression(atop(\x27logit(seq. success):\x27, paste(\x27"$1" coeff. (\x27,alpha^\"" $1 "\", \x27)\x27)))\talpha^"$1}' \
		<(awk -F'\t' '{if ($2=="seq") print $1}' \
			<(grep "\\:" fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_comm_type_design_cols.tsv))) \
	<(awk -F'\t' '{print "logit_prob_MI_coeffs["NR"]\t\"expression(atop(\x27logit(prob. MI):\x27, paste(\x27"$1" coeff. (\x27,beta["NR"], \x27)\x27)))\"\tbeta_"NR}' \
		<(awk -F'\t' '{if ($2=="mi") print $1}' fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_sequencing_technology_comm_type_design_cols.tsv)) |\
	sed '/^$/d' | \
	sed 's/sequencing_technology//g' |	\
	sed 's/bait_capture/bait capture/g' | \
	sed 's/:log10_copies/-vl/g' | \
	sed 's/-vl coeff/ VL coeff/g' | \
	sed 's/alpha^\"bait capture/alpha^\"bait/g' | \
	sed 's/alpha^bait capture/alpha^bait/g' | \
	sed 's/comm_type//g' \
	> config/tmp.tsv

# sort rows
cat \
	<(head -n 1 config/tmp.tsv) \
	<(tail -n +2 config/tmp.tsv | tail -n +2 config/tmp.tsv | grep "logit_prob_seq" | sort) \
	<(grep lambda config/tmp.tsv) \
	<(grep epsilon config/tmp.tsv) \
	<(grep prob_MI config/tmp.tsv | grep -v epsilon | grep -v lambda) \
	> config/sequencing_technology_comm_type_empirical_plot_params.tsv
rm -rf config/tmp.tsv


#### ------------------------------------------------------ ####
#### 4. EXTENDED MODEL WITH PROB. MI DEPENDENT ON COMM_TYPE ####
#### ------------------------------------------------------ ####
Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--inputType full \
	--stan stan/extended_model.stan \
	--miDesignMatrix comm_type \
	--seqDesignMatrix "sequencing_technology" "sequencing_technology:log10_copies" \
	--scaleVars 'log10_copies:sequencing_technology' \
	--designRisks \
		"logit_prob_MI_coeffs[1]+logit_prob_MI" \
		"logit_prob_MI" \
	--designRiskRatios \
		"logit_prob_MI_coeffs[1]+logit_prob_MI:logit_prob_MI" \
	--outAppend comm_type


# need to format plot params file based on design matrix
n_intercept=$(
	awk -F'\t' '{if ($2=="seq") print $0}' fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_comm_type_design_cols.tsv | \
		grep -v "\\:"  | \
		grep -v "col" | wc -l)

cat \
	config/extended_empirical_plot_params_base.tsv \
	<(echo "\n") \
	<(awk -F'\t' '{print "logit_prob_seq_coeffs["NR"]\texpression(atop(\x27logit(seq. success):\x27, paste(\x27"$1" coeff. (\x27,alpha^\"" $1 "\", \x27)\x27)))\talpha^"$1}' \
		<(awk -F'\t' '{if ($2=="seq") print $1}' \
			<(grep -v "\\:" fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_comm_type_design_cols.tsv))) \
	<(awk -F'\t' -v n=$n_intercept '{print "logit_prob_seq_coeffs["NR+n"]\texpression(atop(\x27logit(seq. success):\x27, paste(\x27"$1" coeff. (\x27,alpha^\"" $1 "\", \x27)\x27)))\talpha^"$1}' \
		<(awk -F'\t' '{if ($2=="seq") print $1}' \
			<(grep "\\:" fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_comm_type_design_cols.tsv))) \
	<(awk -F'\t' '{print "logit_prob_MI_coeffs["NR"]\t\"expression(atop(\x27logit(prob. MI):\x27, paste(\x27"$1" coeff. (\x27,beta["NR"], \x27)\x27)))\"\tbeta_"NR}' \
		<(awk -F'\t' '{if ($2=="mi") print $1}' fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_comm_type_design_cols.tsv)) |\
	sed '/^$/d' | \
	sed 's/sequencing_technology//g' |	\
	sed 's/bait_capture/bait capture/g' | \
	sed 's/:log10_copies/-vl/g' | \
	sed 's/-vl coeff/ VL coeff/g' | \
	sed 's/alpha^\"bait capture/alpha^\"bait/g' | \
	sed 's/alpha^bait capture/alpha^bait/g' | \
	sed 's/comm_type//g' \
	> config/tmp.tsv

# sort rows
cat \
	<(head -n 1 config/tmp.tsv) \
	<(tail -n +2 config/tmp.tsv | tail -n +2 config/tmp.tsv | grep "logit_prob_seq" | sort) \
	<(grep lambda config/tmp.tsv) \
	<(grep epsilon config/tmp.tsv) \
	<(grep prob_MI config/tmp.tsv | grep -v epsilon | grep -v lambda) \
	> config/commtype_empirical_plot_params.tsv
rm -rf config/tmp.tsv


#### ----------------------------------------------------------------------------------- ####
#### 5. EXTENDED MODEL WITH PROB. MI DEPENDENT ON COMM_TYPE, SEXPEVER, & SEQUENCING_TECH ####
#### ----------------------------------------------------------------------------------- ####
# ONLY PERFORM THIS ANALYSIS AMONG MEN
Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par_m.tsv \
	--inputType full \
	--stan stan/extended_model.stan \
	--seqDesignMatrix "sequencing_technology" "sequencing_technology:log10_copies" \
	--miDesignMatrix "comm_type" "comm_type:plhiv_sexpever_std" \
	--scaleVars 'log10_copies:sequencing_technology' \
	--miMissingDat \
		'plhiv_sexpever_std:93;3;60;plhiv_sexpeverMany_shape;plhiv_sexpeverMany_scale;plhiv_sexpever_mean' \
		'plhiv_sexpever_std:100;1;60;plhiv_sexpeverMissing_shape;plhiv_sexpeverMissing_scale;plhiv_sexpever_mean' \
	--outAppend sexpever_men

Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par_m_complete.tsv \
	--inputType full \
	--stan stan/extended_model.stan \
	--seqDesignMatrix "sequencing_technology" "sequencing_technology:log10_copies" \
	--miDesignMatrix "comm_type" "comm_type:plhiv_sexpever_std" \
	--scaleVars 'log10_copies:sequencing_technology' \
	--outAppend sexpever_men_complete

# need to format plot params file based on design matrix
n_intercept=$(
	awk -F'\t' '{if ($2=="seq") print $0}' fit/211220_allreads_phsc_all_subgraphs_format_par_m_extended_model_sexpever_men_design_cols.tsv | \
		grep -v "\\:"  | \
		grep -v "col" | wc -l)

cat \
	config/extended_empirical_plot_params_base.tsv \
	<(echo "\n") \
	<(awk -F'\t' '{print "logit_prob_seq_coeffs["NR"]\texpression(atop(\x27logit(seq. success):\x27, paste(\x27"$1" coeff. (\x27,alpha^\"" $1 "\", \x27)\x27)))\talpha^"$1}' \
		<(awk -F'\t' '{if ($2=="seq") print $1}' \
			<(grep -v "\\:" fit/211220_allreads_phsc_all_subgraphs_format_par_m_extended_model_sexpever_men_design_cols.tsv))) \
	<(awk -F'\t' -v n=$n_intercept '{print "logit_prob_seq_coeffs["NR+n"]\texpression(atop(\x27logit(seq. success):\x27, paste(\x27"$1" coeff. (\x27,alpha^\"" $1 "\", \x27)\x27)))\talpha^"$1}' \
		<(awk -F'\t' '{if ($2=="seq") print $1}' \
			<(grep "\\:" fit/211220_allreads_phsc_all_subgraphs_format_par_m_extended_model_sexpever_men_design_cols.tsv))) \
	<(awk -F'\t' '{print "logit_prob_MI_coeffs["NR"]\t\"expression(atop(\x27logit(prob. MI):\x27, paste(\x27"$1" coeff. (\x27,beta["NR"], \x27)\x27)))\"\tbeta_"NR}' \
		<(awk -F'\t' '{if ($2=="mi") print $1}' fit/211220_allreads_phsc_all_subgraphs_format_par_m_extended_model_sexpever_men_design_cols.tsv)) |\
	sed '/^$/d' | \
	sed 's/sequencing_technology//g' |	\
	sed 's/bait_capture/bait capture/g' | \
	sed 's/:log10_copies/-vl/g' | \
	sed 's/-vl coeff/ VL coeff/g' | \
	sed 's/alpha^\"bait capture/alpha^\"bait/g' | \
	sed 's/alpha^bait capture/alpha^bait/g' | \
	sed 's/sexM/men/g' | \
	sed 's/comm_type//g' | \
	sed 's/plhiv_sexpever_std/sexpever/g' \
	> config/tmp.tsv

# sort rows
cat \
	<(head -n 1 config/tmp.tsv) \
	<(tail -n +2 config/tmp.tsv | tail -n +2 config/tmp.tsv | grep "logit_prob_seq" | sort) \
	<(grep lambda config/tmp.tsv) \
	<(grep epsilon config/tmp.tsv) \
	<(grep prob_MI config/tmp.tsv | grep -v epsilon | grep -v lambda) \
	> config/sexpever_men_empirical_plot_params.tsv
rm -rf config/tmp.tsv


#### ----------------------------------------- ####
#### 6. EXTENDED MODEL WITH VARIABLE SELECTION ####
#### ----------------------------------------- ####
Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--inputType full \
	--stan stan/extended_model_hsp.stan \
	--miDesignMatrix "round" "male_circumcision" "comm_type" "sequencing_technology" "comm_type:sex" "comm_type:married" "comm_type:age_cat_coarse" \
		"comm_type:in_migrant" "comm_type:barworker" \
	--seqDesignMatrix "sequencing_technology" "sequencing_technology:log10_copies" \
	--scaleVars 'log10_copies:sequencing_technology' \
	--populationRiskRatios "comm_type='fishing':comm_type='inland'" \
	--outAppend var_select


# need to format plot params file based on design matrix
n_intercept=$(
	awk -F'\t' '{if ($2=="seq") print $0}' fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_hsp_var_select_design_cols.tsv | \
		grep -v "\\:"  | \
		grep -v "col" | wc -l)

cat \
	config/extended_empirical_plot_params_base.tsv \
	<(echo "\n") \
	<(awk -F'\t' '{print "logit_prob_seq_coeffs["NR"]\texpression(atop(\x27logit(seq. success):\x27, paste(\x27"$1" coeff. (\x27,alpha^\"" $1 "\", \x27)\x27)))\talpha^"$1}' \
		<(awk -F'\t' '{if ($2=="seq") print $1}' \
			<(grep -v "\\:" fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_hsp_var_select_design_cols.tsv))) \
	<(awk -F'\t' -v n=$n_intercept '{print "logit_prob_seq_coeffs["NR+n"]\texpression(atop(\x27logit(seq. success):\x27, paste(\x27"$1" coeff. (\x27,alpha^\"" $1 "\", \x27)\x27)))\talpha^"$1}' \
		<(awk -F'\t' '{if ($2=="seq") print $1}' \
			<(grep "\\:" fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_hsp_var_select_design_cols.tsv))) \
	<(awk -F'\t' '{print "w["NR"]\t\"expression(atop(\x27logit(prob. MI):\x27, paste(\x27"$1" coeff. (\x27,beta["NR"], \x27)\x27)))\"\tbeta_"NR}' \
		<(awk -F'\t' '{if ($2=="mi") print $1}' fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_hsp_var_select_design_cols.tsv)) |\
	sed '/^$/d' | \
	sed 's/sequencing_technology//g' |	\
	sed 's/bait_capture/bait capture/g' | \
	sed 's/:log10_copies/-vl/g' | \
	sed 's/-vl coeff/ VL coeff/g' | \
	sed 's/alpha^\"bait capture/alpha^\"bait/g' | \
	sed 's/alpha^bait capture/alpha^bait/g' | \
	sed 's/alpha^bait capture/alpha^bait/g' | \
	sed 's/plhiv_//g' | \
	sed 's/_std//g' | \
	sed 's/TRUE//g' | \
	sed 's/sexM/men/g' | \
	sed 's/_cat_coarse/ /g' | \
	sed 's/male_/male /g' | \
	sed 's/comm_type//g' \
 	> config/tmp.tsv

# sort rows
cat \
	<(head -n 1 config/tmp.tsv) \
	<(tail -n +2 config/tmp.tsv | tail -n +2 config/tmp.tsv | grep "logit_prob_seq" | sort) \
	<(grep lambda config/tmp.tsv) \
	<(grep epsilon config/tmp.tsv) \
	<(grep prob_MI config/tmp.tsv | grep -v epsilon | grep -v lambda) \
	<(grep "w\\[" config/tmp.tsv) \
	> config/var_select_empirical_plot_params.tsv
rm -rf config/tmp.tsv


