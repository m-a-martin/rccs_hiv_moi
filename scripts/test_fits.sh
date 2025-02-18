Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--inputType full \
	--stan stan/deep-phyloMI.stan \
	--filter "id_subgraph_reads > 0 & window_type == 'unique'" \
	--outAppend test

Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--inputType full \
	--stan stan/deep-phyloMI.stan \
	--filter "id_subgraph_reads > 0 & window_type == 'unique'" \
	--outAppend test \
	--seqDesignMatrix "sequencing_technology" "log10_copies" "sequencing_technology*log10_copies"


Rscript scripts/run_stan.R \
	--dat output/211220_allreads_phsc_all_subgraphs_format_par.tsv \
	--inputType full \
	--stan stan/deep-phyloMI.stan \
	--filter "id_subgraph_reads > 0 & window_type == 'unique'" \
	--outAppend test2 \
	--seqDesignMatrix "sequencing_technology" "log10_copies" "sequencing_technology*log10_copies" \
	--miDesignMatrix "comm_type" "sex" "sex*comm_type" \
	--strataPrevs "comm_type == 'fishing'" \
	--multivariateRisks "comm_type == 'fishing'"
	--miCoeffPriors "shrinkage"

	 \
	--miDesignMatrix age_cat_coarse sex comm_type \
	--seqDesignMatrix "sequencing_technology" "sequencing_technology:log10_copies" \
	--scaleVars 'log10_copies:sequencing_technology' \
	--strataRisks "comm_type == 'fishing'" "comm_type == 'inland'" \
	--strataRiskRatios "comm_type == 'fishing':comm_type=='inland'" \
	--outAppend age_sex_comm