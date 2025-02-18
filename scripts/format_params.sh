base=$1
label=$2
data_suffix=$3

cat $base \
	<(echo \n) \
	<(awk -F'\t' '{print "logit_prob_seq_coeffs["NR"]\t\"expression(atop(\x27logit(seq. success):\x27, paste(\x27"$1" coeff. (\x27,alpha["NR"], \x27)\x27)))\"\talpha_"NR}' \
			<(awk -F'\t' '{if ($2=="seq") print $1}' \
				fit/211220_allreads_phsc_all_subgraphs_format_par_${data_suffix}deep-phyloMI_${label}_design_cols.tsv)) \
	<(awk -F'\t' '{print "logit_prob_mi_coeffs["NR"]\t\"expression(atop(\x27logit(prob. MI):\x27, paste(\x27"$1" coeff. (\x27,beta["NR"], \x27)\x27)))\"\tbeta_"NR}' \
			<(awk -F'\t' '{if ($2=="mi") print $1}' \
				fit/211220_allreads_phsc_all_subgraphs_format_par_${data_suffix}deep-phyloMI_${label}_design_cols.tsv)) | \
	sed '/^$/d' | \
	sed 's/sequencing_technology//g' |	\
	sed 's/bait_capture/bait/g' | \
	sed 's/log10_copies/VL/g' | \
	sed 's/age_cat_coarse/age /g' | \
	sed 's/sexM/men/g' | \
	sed 's/sexF/women/g' | \
	sed 's/comm_type//g' | \
	sed 's/plhiv_sexpever_std/sex partners/g' | \
	sed 's/male_circumcisionFALSE/uncircumcised/g' | \
	sed 's/male_circumcisionTRUE/circumcised/g' | \
	sed 's/marriedFALSE/unmarried/g' | \
	sed 's/marriedTRUE/married/g' | \
	sed 's/in_migrantFALSE/non-migrant/g' | \
	sed 's/in_migrantTRUE/migrant/g' | \
	sed 's/barworkerFALSE/non-barworker/g' | \
	sed 's/barworkerTRUE/barworker/g' | \
	sed "s/\*/'%\*%'/g" \
	> config/tmp.tsv

# sort rows
cat \
	<(head -n 1 config/tmp.tsv) \
	<(tail -n +2 config/tmp.tsv | tail -n +2 config/tmp.tsv | grep "logit_prob_seq" | sort) \
	<(grep lambda config/tmp.tsv) \
	<(grep epsilon config/tmp.tsv) \
	<(grep prob_mi config/tmp.tsv | grep -v epsilon | grep -v lambda) \
	> config/empirical_${label}_plot_params.tsv
rm -rf config/tmp.tsv
