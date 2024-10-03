suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(require(patchwork))
suppressMessages(source('scripts/utils.R'))


p <- arg_parser("plot summary of simulated data")
p <- add_argument(p, "--dat", help="input data file", nargs=1)
p = add_argument(p, "--getWindowType1", help="which set of x axis windows to get", default='all', type="character", nargs=1)
p = add_argument(p, "--getWindowType2", help="which set of y axis windows to get", default='unique', type="character", nargs=1)
p <- add_argument(p, "--axesLabels", help="how to label axes", default=c("all", "non-overlapping"), type="character", nargs=2)
p <- add_argument(p, "--out", help="output figure name", nargs=1)
p <- add_argument(p, "--colors", help="tsv file with color codes", nargs=1, default='config/colors.tsv')
args <- parse_args(p)

#args$out="uniq_uniq_alt_corr"
#args$dat="output/211220_allreads_phsc_all_subgraphs_format.tsv"
#args$getWindowType1="unique"
#args$getWindowType2="unique_alt"
#args$axesLabels=c("non-overlapping","non-overlapping alt.")

cols = read_tsv(args$colors, show_col_types=FALSE)
args$colors_dict = setNames(cols$color, cols$var)

#args$dat = 'output/211220_allreads_phsc_all_subgraphs_format_bc.tsv'
if (args$getWindowType1 == "all"){
	x_n_d = summarise_n_d(tabulate_n_d(read_tsv(args$dat, show_col_types = FALSE))) %>%
		select(id, N_obs, MI_obs) %>%
		rename(c('N_obs_x' = 'N_obs', 'MI_obs_x' = 'MI_obs'))
}else{
	x_n_d = summarise_n_d(tabulate_n_d(read_tsv(args$dat, show_col_types = FALSE) %>% 
			filter(window_type == args$getWindowType1))) %>%
		select(id, N_obs, MI_obs) %>%
		rename(c('N_obs_x' = 'N_obs', 'MI_obs_x' = 'MI_obs'))
}

y_n_d = summarise_n_d(tabulate_n_d(read_tsv(args$dat, show_col_types = FALSE) %>% 
		filter(window_type == args$getWindowType2)))  %>%
	select(id, N_obs, MI_obs)  %>%
	rename(c('N_obs_y' = 'N_obs', 'MI_obs_y' = 'MI_obs'))

n_d = x_n_d %>% full_join(y_n_d, by='id') %>%
	mutate(
		N_obs_x = replace_na(N_obs_x, 0),
		MI_obs_x = replace_na(MI_obs_x, 0),
		N_obs_y = replace_na(N_obs_y, 0),
		MI_obs_y = replace_na(MI_obs_y, 0))

p1 = ggplot(n_d, aes(x=N_obs_x, y=N_obs_y)) +
	geom_point(shape=21, color='#333333') +
	ylab(args$axesLabels[2]) +
	xlab(args$axesLabels[1]) +
	ggtitle(expression(paste("sequenced windows (", N[i]^{obs}, ")"))) +
	gtheme


p2 = ggplot(n_d, aes(x=MI_obs_x, y=MI_obs_y)) +
	geom_point(shape=21, color=unname(args$colors_dict["multiple_subgraphs"])) +
	ylab(args$axesLabels[2]) +
	xlab(args$axesLabels[1]) +
	ggtitle(expression(paste("multiple subgraph windows (", MI[i]^{obs}, ")"))) +
	gtheme

p = p1 + p2 + plot_annotation(tag_levels = list(c('a'), '1'))  & 
	theme(plot.tag = element_text(face = 'bold', color='#333333'))

ggsave(paste(c('figures/', args$out,'.pdf'), collapse=''), p, width=12.8, height=4.8)
