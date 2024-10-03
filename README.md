# Identification of HIV multiple infections from deep-sequence within-host phylogenies
This repo contains the code necessary to replicate the data analysis and visualization in the manuscript ''Bayesian identification of HIV multiple infections among people living with HIV in Uganda from whole-genome deep sequence data''. The analysis can be rerun with the `scripts/run_all.sh` wrapper script which in turn calls a number of other wrapper shell scripts and Rscript files. Most of this code was written by Michael A. Martin (post-doctoral research fellow at Johns Hopkins Medicine), who can be reached at [mmart108@jhmi.edu](mailto:mmart108@jhmi.edu) with any questions. 

**Please note that all analysis and text is currently a work in progress**

## Analysis
The steps of the analyis are as follows: 
1. Format raw metadata (`scripts/format_metadata.R`): Formats a number of input metadata files into a consolidated metadata file (`data/input_metadata.tsv`). The input metadata files are not shared publicly due to privacy concerns so this script is shared for reference only. 
2. Process `phyloscanner` output (`scripts/estimate_phsc_moi.R`, `scripts/format_phsc_moi.R`): Takes as input a number of `.rda` output files from `phyloscanner` and processes them to identify the number of subgraphs in each window for each processed sample and the cophenetic distance between subgraph MRCAs among sample windows with multiple subgraphs. `Phyloscanner` output is then combined with metadata to generate a single consolidated data file, `output/211220_allreads_phsc_all_subgraphs_format.tsv`. Similarly to the above, phyloscanner outputs are not shared publicly but these scripts are shared for reference. 
3. Get representative phylogenetic trees (`scripts/get_tree_dat`): Gets representative phylogenetic tree objects from `phyloscanner` output files for empirical data summary figure. 
4. Simulation study (`scripts/simulation_study.sh`): Simulates a variety of data sets using the `scripts/simulation_model.R` file and fits our inference models (`stan/*.stan`) using `cmdstanR` (`scripts/run_stan.R`). All model fits are saved in `fit/`. Fits are not stored on GitHub given file size restraints but will be uploaded to Zenodo prior to submission. 
5. Empirical model fits (`scripts/empirical_model_fits.sh`): Runs a variety of inference models on our empirical data file generated in Step #2. All model fits are saved in `fit/`. 
6. Plot figures (`scripts/plot_empirical_data.sh`, `scripts/plot_simulation_model_fits.sh`, `scripts/scripts/plot_empirical_model_fits.sh`): First plots a summary of our empirical data, then generates figures from our simulation fits, and finally from our model fits to empirical data. All figures are saved in `figures/` in PDF format. 
7. Prepare supplementary files: Simple bash script to process some of the `phyloscanner` output for File S1. 
8. Calculate statistics (`scripts/calc_stats.R`): Takes as input a number of data and model fit files and generates a file with statistics (`output/statistics.csv`) used in rendering the manuscript. 
9. Render document with `lualatex` and `bibtex`.  

## Dependencies
The workflow above requires a number of dependencies. The program versions provided are those which were used. Others may work. 
- adephylo v.1.1-16
- ape v.5.8
- argparser v.0.7.2
- bayesplot v.1.11.1
- bibtex v.0.99d
- conflicted v.1.2.0
- cowplot v.1.1.3
- fitdistrplus v.1.1-11
- ggplot2 v.3.5.1
- ggtree v.3.12.0
- haven v.2.5.4
- HDInterval v.0.2.4
- luahbtex v.1.15.0
- patchwork v.1.2.0
- phytools v.2.1-1
- posterior v.1.6.0
- R v.4.4.1
- stringr v.1.5.1
- tidytree v.0.4.6
- tidyverse (dplyr v.1.1.4, forcats v.1.0.0, lubridate v.1.9.3, purr v.1.0.2, readr v.2.1.5, tibble v.3.2.1, tidyr v.1.3.1)




