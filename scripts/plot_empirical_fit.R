suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(require(cowplot))
suppressMessages(require(patchwork))
suppressMessages(require(HDInterval))
suppressMessages(source('scripts/utils.R'))

p <- arg_parser("plot summary of simulated data")
p <- add_argument(p, "--dat", help="input data file", nargs=1)
p <- add_argument(p, "--fit", help="stan fit object", nargs=1)
p <- add_argument(p, "--out", help="output figure name", nargs=1)
p <- add_argument(p, "--colors", help="tsv file with color codes", nargs=1, default='config/colors.tsv')
args <- parse_args(p)

#args$dat = 'output/211220_allreads_phsc_all_subgraphs_format_par.tsv'
#args$fit = 'fit/211220_allreads_phsc_all_subgraphs_format_par_full_model_.Rds'
#args$colors = 'config/colors.tsv'

cols = read_tsv(args$colors, show_col_types=FALSE)
args$colors_dict = setNames(cols$color, cols$var)

dat = summarise_n_d(tabulate_n_d(read_tsv(args$dat, show_col_types = FALSE) %>% 
	filter(window_type == 'unique'))) 
fit = readRDS(args$fit)

p1 = plot_ind_posterior_pred(dat,fit, args$colors_dict, guide=TRUE) + xlab('participants')

param='prob_MI'
p_dat = as_tibble(fit$draws(
		variables = c(param),
		  inc_warmup = FALSE,
		  format = "draws_df")) %>%
	select(all_of(param))
median_prob_mi = median(p_dat$prob_MI)
# get HPD bounds
bounds95 = hdi(p_dat, credMass=0.95)[,1]
bounds50 = hdi(p_dat, credMass=0.50)[,1]
bounds = c(-Inf, bounds95[1], bounds50[1], bounds50[2], bounds95[2], Inf)
p2 = plot_shadded_hist(p_dat, bounds, approx_bins=50,
			cols=c(unname(args$colors_dict["mi_hpd3"]), 
				unname(args$colors_dict["mi_hpd2"]), 
				unname(args$colors_dict["mi_hpd1"]))) +
	xlab(expression(atop('prevalence of', paste('multiple infections (', bar(delta), ')'))))

n_mi_cutoff_dat = read_tsv(str_replace(args$fit, '.Rds', '_n_mi_dichotomized.tsv'), show_col_types=FALSE)

n_mi_cutoff = n_mi_cutoff_dat %>% filter(cutoff != 'empirical' & cutoff != 'n_mi') %>%
	mutate(cutoff = as.numeric(cutoff))

p3 = ggplot(n_mi_cutoff, aes(x=cutoff, y=`0.50`)) +
	geom_ribbon(aes(ymin=`0.025`, ymax=`0.975`), fill=args$colors_dict["mi_hpd3"]) +
	geom_ribbon(aes(ymin=`0.25`, ymax=`0.75`), fill=args$colors_dict["mi_hpd2"]) +
	geom_line(color=args$colors_dict["mi_hpd1"]) +
	geom_hline(aes(yintercept=(n_mi_cutoff_dat %>% filter(cutoff == 'n_mi'))$`0.50`), 
		linetype='dashed', color='#333333') +
	ylab('participants w/\nmultiple infection') +
	xlab('posterior prob. of\nmultiple infection threshold') +
	ylim(0,NA) +
	xlim(0,NA) +
	gtheme

#p = p1 +  plot_spacer() + p2 + 
#	plot_layout(widths=c(1, 0.025, 0.6)) + 
#	plot_annotation(tag_levels = list(c('a'), '1'))  & 
#	theme(plot.tag = element_text(face = 'bold', color='#333333'))

p = p1 / (p3 + plot_spacer() + p2 + plot_layout(widths=c(1, 0.05, 1))) +
	plot_annotation(tag_levels=list(c('a'), '1')) & 
	theme(plot.tag = element_text(face = 'bold', color='#333333'))

ggsave(paste(c('figures/', args$out,'.pdf'), collapse=''), p, width=10, height=8)


