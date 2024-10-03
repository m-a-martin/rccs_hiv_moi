#!/usr/bin/env bash

#### ---------------------- ####
#### 1. FORMAT RAW METADATA ####
#### ---------------------- ####
# this step is internal and processes data
# that will not be shared publicly 
Rscript scripts/format_metadata.R


#### ------------------------------ ####
#### 2. PROCESS PHYLOSCANNER OUTPUT ####
#### ------------------------------ ####
mkdir output
mkdir logs
# loop through rda objects, get # of subgraphs
# takes an hour or two to run, can do it in parallel
# and concatenate output to speed things up
Rscript scripts/estimate_phsc_moi.R \
		output/211220_allreads_phsc_all_subgraphs.tsv \
		data/211220_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd/ptyr*.rda \
		> logs/logs_allreads_phsc.out \
		2>&1

# format output, add in metadata 
Rscript scripts/format_phsc_moi.R \
	--dat 'output/211220_allreads_phsc_all_subgraphs.tsv' \
	--ids 'output/211220_allreads_phsc_all_subgraphs_ids.tsv' \
	--metadata 'data/input_metadata_internal.tsv' 
# some additional filtering of input data
bash scripts/empirical_datasets.sh


#### ------------------------------------ ####
#### 3. REPRESENTATIVE PHYLOGENETIC TREES ####
#### ------------------------------------ ####
# used in empirical data summary figure
Rscript scripts/get_tree_dat.R \
	--dat data/211220_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd/ptyr298_workspace.Rda \
	--windowStart 1800 \
	--sample AID2642-fq2

Rscript scripts/get_tree_dat.R \
	--dat data/211220_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd/ptyr12_workspace.Rda \
	--windowStart 1800 \
	--sample AID2211-fq1


# PUBLICLY AVAILABLE ANALYSIS BEGINS HERE
#### ------------------- ####
#### 4. SIMULATION STUDY ####
#### ------------------- ####
bash scripts/simulation_study.sh


#### ------------------------------ ####
#### 5. FIT MODEL TO EMPIRICAL DATA ####
#### ------------------------------ ####
bash scripts/empirical_model_fits.sh
# sensitivity to choice of gneome windows
bash scripts/genome_window_sensitivity.sh


#### --------------- ####
#### 6. PLOT FIGURES ####
#### --------------- ####
# first empirical data
bash scripts/plot_empirical_data.sh
# then simulation fits
bash scripts/plot_simulation_model_fits.sh
# then empirical fits
bash scripts/plot_empirical_model_fits.sh


#### ------------------------------ ####
#### 7. PREPARE SUPPLEMENTARY FILES ####
#### ------------------------------ ####
awk -F'\t' '{print $2}' output/211220_allreads_phsc_all_subgraphs_refs.tsv | \
	sort | uniq | sed 's/REF_//g' | grep -v ref | grep -v "^$" \
	> manuscript/supplementary_files/file_s1.txt


#### --------------------------------------- ####
#### 8. CALCUALTE STATISTICS FOR TEXT/TABLES ####
#### --------------------------------------- ####
Rscript scripts/calc_stats.R \
	--phscDat 'output/211220_allreads_phsc_all_subgraphs_format_par.tsv' \
	--metadata 'output/211220_allreads_phsc_metadata.tsv' \
	--baseBaseFit 'fit/base_simulation_base_model__summary.tsv' \
	--baseParams 'simulations/base_simulation_params.tsv' \
	--fullBaseFit 'fit/full_simulation_base_model__summary.tsv' \
	--fullFullFit 'fit/full_simulation_full_model__summary.tsv' \
	--extParams 'simulations/extended_simulation_params.tsv' \
	--extExtFit 'fit/extended_simulation_extended_model_hsp__summary.tsv' \
	--fullParams 'simulations/full_simulation_params.tsv' \
	--empiricalFullFit 'fit/211220_allreads_phsc_all_subgraphs_format_par_full_model__summary.tsv' \
	--empiricalFullAltFit 'fit/211220_allreads_phsc_all_subgraphs_format_par_full_model_alt_summary.tsv' \
	--empiricalSeqFit 'fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_sequencing_technology_summary.tsv' \
	--empiricalCommFit 'fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_comm_type_summary.tsv' \
	--empiricalSexpeverMenFit 'fit/211220_allreads_phsc_all_subgraphs_format_par_m_extended_model_sexpever_men_summary.tsv' \
	--empiricalSexpeverMenCompleteFit 'fit/211220_allreads_phsc_all_subgraphs_format_par_m_complete_extended_model_sexpever_men_complete_summary.tsv' \
	--empiricalVarSelectFit 'fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_hsp_var_select_summary.tsv'

# format database file, leads to faster latex compilation
/Applications/datatooltk/bin/datatooltk --csv output/statistics.csv --output output/statistics.dbtex

#### --------------------- ####
#### 9. COMPILE MANUSCRIPT ####
#### --------------------- ####
# copy over from overleaf directory 
cp ~/Library/CloudStorage/Dropbox/Apps/Overleaf/HIV\ multiple\ infection\ manuscript/manuscript/*.tex \
	manuscript/.
cp ~/Library/CloudStorage/Dropbox/Apps/Overleaf/HIV\ multiple\ infection\ manuscript/manuscript/*.bst \
	manuscript/.
cp ~/Library/CloudStorage/Dropbox/Apps/Overleaf/HIV\ multiple\ infection\ manuscript/manuscript/*.bib\
	manuscript/.
cp ~/Library/CloudStorage/Dropbox/Apps/Overleaf/HIV\ multiple\ infection\ manuscript/manuscript/supplementary_files/*.tex \
	manuscript/supplementary_files/.

# todo equation numbering
cd manuscript
rm -rf ./*_vals_figs.tex
rm -rf ./*_vals.tex
rm -rf ./tmp.tex
# create a version with no figures
for i in ./*.tex; do
	bash ../scripts/format_manuscript.sh $i ../output/statistics.scsv
done

# do the same for the supplement
cd supplementary_files
rm -rf ./*_vals_figs.tex
rm -rf ./*_vals.tex
rm -rf ./tmp.tex
scp ../references.bib .
scp ../plos2015.bst .
# create a version with no figures
for i in ./*.tex; do
	bash ../../scripts/format_manuscript.sh $i ../../output/statistics.scsv
done

cd ..

# move to sharing folder
mkdir -p to_share
mkdir -p to_share/word_figs
mkdir -p to_share/word_nofigs
mkdir -p to_share/pdf

scp *.pdf to_share/pdf/.
scp supplementary_files/*.pdf to_share/pdf/.
scp supplementary_files/*.txt to_share/pdf/.
scp supplementary_files/*.txt to_share/word_figs/.
scp supplementary_files/*.txt to_share/word_nofigs/.
scp supplementary_files/*.csv to_share/pdf/.
scp supplementary_files/*.csv to_share/word_figs/.
scp supplementary_files/*.csv to_share/word_nofigs/.
scp *.docx to_share/word_nofigs/.
scp supplementary_files/*.docx to_share/word_nofigs/.
for i in to_share/word_nofigs/*_figs.docx; do
	mv $i ${i/nofigs/figs}
done

cd ..
#### end final code here ####
