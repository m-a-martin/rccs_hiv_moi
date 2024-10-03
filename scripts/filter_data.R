suppressMessages(library(tidyverse))
suppressMessages(require(readxl))

exclude = c('RK-D119257', 'RK-F125812', 'RK-B113255')

hiv_dr = tibble()
for (sheet in c("5%_10reads", "5%_3reads", "20%_10reads", "20%_3reads")){
	hiv_dr = bind_rows(
		hiv_dr, 
		read_excel("data/Rakai_Drug_Resistance_20220921.xlsx", sheet=sheet)  %>%
			mutate(
				sampleID_clean = str_split(sampleID, '5pct', simplify=TRUE)[,1],
				sheet=sheet))
}

# read in mapping to rakai ID and date
metadata = read_csv('data/rakai_sequence_id_mappings_cohort.csv', show_col_types=FALSE) %>%
	mutate(drmseq_prefix = str_split(drmseq_prefix, '20pct', simplify=TRUE)[,1])

exclude_seqs = (metadata %>% filter(pt_id %in% exclude))$drmseq_prefix

hiv_dr_filter = hiv_dr %>% filter(!(sampleID_clean %in% exclude_seqs))

write_tsv(hiv_dr_filter %>% select(-sampleID_clean), "data/Rakai_Drug_Resistance_20220921_filter.tsv")