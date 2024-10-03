#!/usr/bin/env bash

# filter data for just men 
cat \
	<(head -n 1 output/211220_allreads_phsc_all_subgraphs_format_par.tsv) \
	<(awk -F'\t' 'FNR==NR{split($0, a, "\t"); for (i in a) cols[a[i]] = i; next}{split($0,b,"\t"); if (b[cols["sex"]] == "M") print $0}' \
		<(head -n 1 output/211220_allreads_phsc_all_subgraphs_format_par.tsv) \
		output/211220_allreads_phsc_all_subgraphs_format_par.tsv) \
	> output/211220_allreads_phsc_all_subgraphs_format_par_m.tsv

# check if there are there any sexpever == 92 responses?
awk -F'\t' 'FNR==NR{split($0, a, "\t"); for (i in a) cols[a[i]] = i; next}{split($0,b,"\t"); if (b[cols["plhiv_sexpever_std"]] == 92) print $0}' \
		<(head -n 1 output/211220_allreads_phsc_all_subgraphs_format_par_m.tsv) \
		output/211220_allreads_phsc_all_subgraphs_format_par_m.tsv

# complete case analysis
# filter data for just men with numeric reseponse 
awk -F'\t' 'FNR==NR{split($0, a, "\t"); for (i in a) cols[a[i]] = i; next}{split($0,b,"\t"); if (b[cols["plhiv_sexpever_std"]] != "100" && b[cols["plhiv_sexpever_std"]] != "93") print $0}' \
	<(head -n 1 output/211220_allreads_phsc_all_subgraphs_format_par.tsv) \
	output/211220_allreads_phsc_all_subgraphs_format_par_m.tsv \
	> output/211220_allreads_phsc_all_subgraphs_format_par_m_complete.tsv
