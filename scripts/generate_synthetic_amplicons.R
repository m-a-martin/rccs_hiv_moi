suppressMessages(library(tidyverse))
suppressMessages(library(argparser))
suppressMessages(source('scripts/utils.R'))


p <- arg_parser("generate synthetic amplicon data")
p <- add_argument(p, "--dat", help="input data file", nargs=1)
p <- add_argument(p, "--starts", help="start location of each region", nargs=Inf)
p <- add_argument(p, "--ends", help="end location of each region", nargs=Inf)
p <- add_argument(p, "--filter", default="id_subgraph_reads > 0", help='filter data')
args <- parse_args(p)

#args$dat = 'output/211220_allreads_phsc_all_subgraphs_format_par.tsv'
#args$starts = c(1427, 7941)
#args$ends = c(1816, 8264)

args$starts = as.integer(args$starts)
args$ends = as.integer(args$ends)


dat = tabulate_n_d(read_tsv(args$dat, show_col_types=FALSE) %>% filter(eval(parse(text=args$filter))))
all_windows = dat %>% select(window_start, window_end) %>% unique()
get_windows = tibble()
for (i in 1:length(args$starts)){
	start = args$starts[i]
	end = args$ends[i]
	get_windows = bind_rows(
		get_windows, 
		all_windows %>% filter(window_start >= start & window_end < end) %>%
			arrange(window_start) %>%
			filter(row_number() == 1 | row_number() == n()) %>%
			mutate(amplicon_idx = i))
}

dat = dat %>% left_join(get_windows, by=c('window_start'))

amplicon_dat = dat %>% filter(!is.na(amplicon_idx)) %>%
	group_by(study_id, round, sex, age_cat_coarse, comm_type, amplicon_idx, sequencing_technology, log10_copies, log10_copies_shared) %>%
	summarise(n = 1*(sum(n)==2), d = 1*(sum(d*n) > 0)) %>%
	rename(window_start = amplicon_idx) %>%
	mutate(n_subgraphs = NA, window_end = NA, subgraph_reads_all = NA, subgraph_freq_all = NA, id_subgraph_reads = NA,
		id_subgraph_reads_all = NA,  sg=NA, id_min_sg_dist = NA, min_sg_dist = NA,
		window_type='unique', N_windows=2)

write_tsv(amplicon_dat, 'output/211220_allreads_phsc_all_subgraphs_format_par_amplicon.tsv')
