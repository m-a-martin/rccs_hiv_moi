suppressMessages(library(tidyverse))
suppressMessages(source('scripts/utils.R'))

# MI data
full_d = read_tsv('output/211220_allreads_phsc_all_subgraphs.tsv', show_col_types=FALSE) 

# input data
input = read_tsv('data/210120_RCCSUVRI_phscinput_samples.tsv', show_col_types=FALSE)

# exclusion data
exclusion = c('RK-D119257', 'RK-F125812', 'RK-B113255')
exclusion_id = (input %>% filter(UNIT_ID %in% exclusion))$RENAME_ID

clean_full_d = full_d %>% filter(!(id %in% exclusion_id))

write_tsv(clean_full_d, 'output/211220_allreads_phsc_all_subgraphs.tsv')

# list of IDs
full_id = read_tsv('output/211220_allreads_phsc_all_subgraphs_ids.tsv', show_col_types=FALSE) 
clean_full_id = full_id %>% filter(!(id %in% exclusion_id))
write_tsv(clean_full_id, 'output/211220_allreads_phsc_all_subgraphs_ids.tsv')
