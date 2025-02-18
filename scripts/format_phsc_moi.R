suppressMessages(library(tidyverse))
suppressMessages(library(argparser))
suppressMessages(source('scripts/utils.R'))


p <- arg_parser("format phsc output")
p <- add_argument(p, "--dat", help="input data file", nargs=1)
p <- add_argument(p, "--ids", help="all ids in phsc runs", nargs=1)
p <- add_argument(p, "--metadata", help="metadata", nargs=1)
args <- parse_args(p)

viremic = 3
#in_args = c('output/211220_allreads_phsc_all_subgraphs.tsv', 'output/211220_allreads_phsc_all_subgraphs_ids.tsv', 'data/formatted_metadata.tsv')
#args$dat = 'output/211220_allreads_phsc_all_subgraphs.tsv'
#args$ids = 'output/211220_allreads_phsc_all_subgraphs_ids.tsv'
#args$metadata = 'data/input_metadata_internal.tsv'

in_split = str_split(args$dat, '_', simplify=TRUE)

full_d = read_tsv(args$dat, show_col_types=FALSE) 

# for each sequencing run, get the phsc run with the 
# most subgraphs and most dual subgraphs
keep_runs = full_d %>% group_by(id, window_start, window_end, run) %>%
			summarise(n_subgraphs = length(unique(split.ids)), .groups='drop') %>%
			group_by(run, id) %>%
			summarise(n=n(), d=sum(n_subgraphs > 1), .groups="drop") %>%
			arrange(-n, -d) %>%
			group_by(id) %>%
			slice(1) %>%
			select(run, id)

keep_runs = setNames(keep_runs$run, keep_runs$id)

d = full_d %>% filter(
	keep_runs[full_d$id] == run)

# add in metadata
metadata = read_tsv(args$metadata, show_col_types=FALSE)

# merge metadata and subgraph data
# filter for just viremic samples	
d = metadata %>% left_join(d, by='id') %>%
	filter(!is.na(log10_copies) &
		log10_copies >= viremic & 
		!is.na(run))

# rename split.ids columns
d = d %>% rename(sg = split.ids) %>%
	# add "subject" columns
	mutate(subject = str_split(id, "-", simplify=T)[,1])
	
# label windows if they are unique
# get just non-overlapping windows
window_size=(d$window_end[1] - d$window_start[1]) + 1
windows = seq(min(d$window_start), max(d$window_start) + window_size,window_size)
# remove windows that were filtered out
windows = windows[windows %in% d$window_start]
# get alternate windows
windows_alt = seq(min(d$window_start) + window_size/2, max(d$window_start) + window_size,window_size)
# remove windows that were filtered out
windows_alt = windows_alt[windows_alt %in% d$window_start]

d = d %>% mutate(window_type = case_when(
	window_start %in% windows ~ 'unique', 
	window_start %in% windows_alt ~ 'unique_alt',
	!(window_start %in% windows) & !(window_start %in% windows_alt) ~ 'overlapping'),
	N_windows = length(windows))


# filter visits for just the highest N_obs among unique windows
# if ties, get maximum viral load
par_d = d %>% inner_join(
	summarise_n_d(tabulate_n_d(d %>% filter(window_type == 'unique' & id_subgraph_reads > 0))) %>% 
		select(study_id, round, N_obs, log10_copies) %>%
		group_by(study_id) %>%
		filter(N_obs == max(N_obs)) %>%
		filter(log10_copies == max(log10_copies)),
	by=c('study_id', 'round', 'log10_copies'))

# read in list of ids
ids = read_tsv(args$ids, show_col_types=FALSE) %>%
	unique() %>% 
	mutate(phsc = TRUE)

metadata = 
	metadata %>% 
		left_join(
			ids %>% 
				select(id, phsc) %>% 
				unique(), 
			by=c('id' = 'id')) %>%
		left_join(par_d %>% select(study_id, round, run) %>% unique() %>%
			mutate(phsc_par = TRUE)) %>%
	mutate(phsc = replace_na(phsc, FALSE),
		phsc_par = replace_na(phsc_par, FALSE))


# add phsc dates
phsc_dates = metadata %>% 
	filter(log10_copies >= 3 & phsc == TRUE) %>%
	summarise(
		min_phsc_date = 
			paste(c(
				format(min(visit_dt), "%B"), 
				" ",
				format(min(visit_dt), "%Y")), collapse=''),
		max_phsc_date = 
			paste(c(
				format(max(visit_dt), "%B"), 
				" ",
				format(max(visit_dt), "%Y")), collapse=''))

metadata$min_phsc_date = phsc_dates$min_phsc_date
metadata$max_phsc_date = phsc_dates$max_phsc_date


write_tsv(metadata %>% select(-rccs_study_id, -visit_dt),
	paste(c(in_split[1:3], 'metadata.tsv'), collapse='_'))
write_tsv(d %>% select(-rccs_study_id, -visit_dt),
	paste(c(in_split[1:3], 'all_subgraphs_format.tsv'), collapse='_'))
write_tsv(par_d  %>% select(-rccs_study_id, -visit_dt), 
	paste(c(in_split[1:3], 'all_subgraphs_format_par.tsv'), collapse='_'))
