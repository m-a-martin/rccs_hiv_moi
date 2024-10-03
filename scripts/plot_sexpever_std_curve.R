suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(source('scripts/utils.R'))


p <- arg_parser("plot summary of empirical data")
p <- add_argument(p, "--dat", help="input data file", nargs=1)
p <- add_argument(p, "--colors", help="tsv file with color codes", default='config/colors.tsv', nargs=1)
p <- add_argument(p, "--out", help="output figure name", nargs=1)
args <- parse_args(p)
cols = read_tsv(args$colors, show_col_types=FALSE)
args$colors_dict = setNames(cols$color, cols$var)


#args$dat = 'data/input_metadata.tsv'

# drop ambiguous responses (92 and 93)
# drop 0s
dat = read_tsv(args$dat, show_col_types = FALSE) %>% 
	filter(sex == 'M' & sexpever != 92 & sexpever != 93 & sexpever != 0 & finalhiv == 'P') %>%
	mutate(sex = if_else(sex == 'M', 'men', 'women'),
		comm_type = ordered(comm_type, levels=c('inland', 'fishing')))

mean_est_sexpever = dat %>% filter(finalhiv == 'P') %>%
	select(sex, comm_type, age_cat_fine, plhiv_sexpever_mean) %>% unique()


p = ggplot(mean_est_sexpever, aes(x=age_cat_fine, y=plhiv_sexpever_mean)) +
	geom_bar(stat="identity", fill='#eaeaea', color='#333333') +
	facet_grid(cols = vars(comm_type)) +
	ylab(expression(est.~mean~lifetime~sex~partners~(bar(P[s])))) +
	xlab('age') + 
	gtheme



ggsave(paste(c('figures/', args$out, '.pdf'), collapse=''), p, width=12, height=4.8)



