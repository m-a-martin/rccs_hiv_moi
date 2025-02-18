suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(require(patchwork))
suppressMessages(source('scripts/utils.R'))


p <- arg_parser("plot summary of simulated data")
p <- add_argument(p, "--dat", help="input data file", nargs=1)
p = add_argument(p, "--filter", default="TRUE", help="input data filter", nargs=1)
p <- add_argument(p, "--out", help="output figure name", nargs=1)
p <- add_argument(p, "--colors", help="tsv file with color codes", nargs=1, default='config/colors.tsv')
args <- parse_args(p)
#args$dat = 'output/211220_allreads_phsc_all_subgraphs_format.tsv'
cols = read_tsv(args$colors, show_col_types=FALSE)
args$colors_dict = setNames(cols$color, cols$var)


dat = tabulate_n_d(read_tsv(args$dat, show_col_types = FALSE) %>% 
	filter(eval(parse(text=args$filter)))) %>%
	mutate(vl_cat = cut(log10_copies, breaks=c(3,3.5,4,4.5,5,Inf)))



per_sample = summarise_n_d(dat) %>% 
	group_by(vl_cat, N_obs, sequencing_technology) %>%
	summarise(n_samples=n(), .groups='drop')

idx_var = colnames(per_sample)[
	which(colnames(per_sample) != 'N_obs' & 
		colnames(per_sample) != 'n_samples')]

per_sample = per_sample %>% 
	pivot_wider(names_from=N_obs, values_from=n_samples, values_fill=0) %>%
	pivot_longer(-all_of(idx_var), names_to='N_obs', values_to='n_samples') %>%
	mutate(N_obs = as.numeric(N_obs))

per_window = dat %>% 
	group_by(window_start, window_end, vl_cat, sequencing_technology) %>%
	summarise(n_sequenced = sum(n_subgraphs > 0), n=n(), .groups='drop') %>%
	mutate(window_mid = ( ( window_end + window_start ) / 2 ) ) %>%
	arrange(vl_cat, window_mid)

p1 = ggplot(per_sample, aes(x=N_obs, y=n_samples,
		group=sequencing_technology, fill=sequencing_technology)) + 
	geom_bar(stat="identity", position="dodge", color='#333333') + 
	facet_wrap(~vl_cat, ncol=1) +
	xlab(expression(paste("sequenced windows (", N[i]^{obs}, ")"))) + 
	ylab('number of samples') +
	scale_fill_manual(values=args$colors_dict, name=NULL) +
	gtheme


p2 = ggplot(per_window, aes(x=window_mid, y=n_sequenced, 
  		group=sequencing_technology, fill=sequencing_technology)) + 
	geom_bar(stat="identity", position="dodge", color='#333333', ) +
	facet_wrap(~vl_cat, ncol=1) +
	xlab('position (nt)') + 
	ylab('sequenced samples') +
	xlim(0,NA) +
	scale_fill_manual(values=args$colors_dict, name=NULL) +
	gtheme

p = p1 + p2 + plot_layout(guides = "collect") + plot_annotation(tag_levels = list(c('a'), '1'))  & 
	theme(legend.position = "bottom", plot.tag = element_text(face = 'bold', color='#333333')) 

ggsave(paste(c('figures/', args$out,'.pdf'), collapse=''), p, width=8, height=10)



