suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(require(tidytree))
suppressMessages(source('scripts/utils.R'))


p <- arg_parser("plot summary of empirical data")
p <- add_argument(p, "--dat", help="input data file", nargs=1)
p <- add_argument(p, "--out", help="output figure name", nargs=1)
p <- add_argument(p, "--colors", help="tsv file with color codes", default='config/colors.tsv', nargs=1)
args <- parse_args(p)

#args$dat = "output/211220_allreads_phsc_all_subgraphs_format.tsv"
#args$out = 'test'

dat = read_tsv(args$dat, show_col_types=FALSE)


mrca_bl_dat = dat %>% select(study_id, round, window_start, mrca_dist) %>% unique()

p = ggplot(mrca_bl_dat, aes(x=mrca_dist)) +
	geom_histogram(binwidth=0.005, color='#333333', fill='#eaeaea') +
	scale_x_continuous(limits=c(-0.0025, 0.7), expand=expansion(0.01)) +
	xlab('subs/site') +
	scale_y_log10(breaks=NULL, name='density', expand = expansion(mult = c(0, .15))) +
	gtheme

ggsave(paste(c('figures/', args$out, '.pdf'), collapse=''), p, width=8, height=4.8)
