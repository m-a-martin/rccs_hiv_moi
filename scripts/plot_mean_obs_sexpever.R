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

dat = read_tsv(args$dat, show_col_types = FALSE) %>% 
	filter(
		sexpever != 93 &
		sexpever != 97 & 
		sexpever != 98 & 
		sexpever != 99) %>%
	mutate(hiv=
			ordered(
				if_else(finalhiv == 'P', 'living w/ HIV', 'living w/o HIV'),
				levels=c('living w/o HIV', 'living w/ HIV')),
		comm_type = ordered(comm_type, levels=c('inland', 'fishing')),
		sex = if_else(sex == 'F', 'women', 'men'))

sum_dat = dat %>%
	group_by(hiv, sex, age_cat_fine, comm_type) %>%
	summarise(mean_sexpever = mean(sexpever), .groups='drop')

p = ggplot(sum_dat, aes(x=age_cat_fine, y=mean_sexpever, group=hiv, fill=hiv))+
	facet_grid(cols=vars(comm_type), rows=vars(sex)) + 
	geom_bar(stat="identity", position="dodge", color='#333333')+
	scale_fill_manual(values=args$colors_dict, name=NULL) +
	guides(fill = guide_legend(position = "inside")) +
	xlab('age') +
	ylab('lifetime sex partners') +
	scale_y_continuous(expand = expansion(mult = c(0, .1))) + 
	gtheme +
	theme(legend.justification.inside=c(0,1),
		axis.text.x = element_text(size=10))

ggsave(paste(c('figures/', args$out, '.pdf'), collapse=''), p, width=8, height=8)

