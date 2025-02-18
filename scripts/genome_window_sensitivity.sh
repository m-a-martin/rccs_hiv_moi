#!/usr/bin/env bash

#### ------------------------------------------------------------------ ####
#### 1.EVALUATE DATA CORRELATION BETWEEN ALL WINDOWS AND UNIQUE WINDOWS ####
#### ------------------------------------------------------------------ ####
Rscript scripts/plot_all_uniq_corr.R \
	--out all_uniq_corr \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--filter "id_subgraph_reads > 0"


#### --------------------------------------------------------------------------------- ####
#### 2.EVALUATE DATA CORRELATION BETWEEN UNIQUE WINDOWS AND ALTERNATIVE UNIQUE WINDOWS ####
#### --------------------------------------------------------------------------------- ####
Rscript scripts/plot_all_uniq_corr.R \
	--out uniq_uniq_alt_corr \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--filter "id_subgraph_reads > 0" \
	--getWindowType1 "unique" \
	--getWindowType2 "unique_alt" \
	--axesLabels "non-overlapping" "non-overlapping alt."


#### ---------------------------------------------- ####
#### 3.FIT FULL MODEL TO ALTERNATIVE UNIQUE WINDOWS ####
#### ---------------------------------------------- ####
Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--inputType full \
	--stan stan/deep-phyloMI.stan \
	--miDesignMatrix age_cat_coarse sex comm_type \
	--seqDesignMatrix "sequencing_technology" "log10_copies_shared" "sequencing_technology*log10_copies" \
	--scaleVars 'log10_copies_shared' 'log10_copies:sequencing_technology' \
	--filter "window_type == 'unique_alt' & id_subgraph_reads > 0" \
	--outAppend age_sex_comm_alt


