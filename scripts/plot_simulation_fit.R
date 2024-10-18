suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(require(cowplot))
suppressMessages(require(HDInterval))
suppressMessages(source('scripts/utils.R'))




p <- arg_parser("plot summary of simulated data")
p <- add_argument(p, "--dat", help="input data file", nargs=1)
p <- add_argument(p, "--fit", help="stan fit object", nargs=1)
p <- add_argument(p, "--params", help="simulation parameters", nargs=1)
p <- add_argument(p, "--plotParams", help="parameters to plot", nargs=1)
p <- add_argument(p, "--ncolParams", help="parameters to plot", nargs=1, default=3)
p <- add_argument(p, "--out", help="output figure name", nargs=1)
p <- add_argument(p, "--colors", help="tsv file with color codes", nargs=1, default='config/colors.tsv')
args <- parse_args(p)
#args$fit = 'fit/full_simulation_full_model_.Rds'
#args$dat = 'simulations/full_simulation.tsv'
#args$params = 'simulations/full_simulation_params.tsv'
#args$plotParams = 'config/tmp.tsv'
#args$ncolParams = 5

cols = read_tsv(args$colors, show_col_types=FALSE)
args$colors_dict = setNames(cols$color, cols$var)

plot_params = read_tsv(args$plotParam, show_col_types=FALSE)
dat = read_tsv(args$dat, show_col_types = FALSE) 
fit = readRDS(args$fit)

params = read_tsv(args$params, show_col_types=FALSE) %>%
	filter(arg != 'out')
params_dict = setNames(as.numeric(params$value), params$arg)

p1 = plot_n_d_scatter(dat, args$colors_dict, grouping="truth")
p2 = plot_ind_posterior_pred(dat, fit, args$colors_dict, grouping="truth")

plot_list = list()
for (idx in seq(1,nrow(plot_params))){
	param = plot_params$param[idx]
	p_dat = as_tibble(fit$draws(
			variables = c(param),
			  inc_warmup = FALSE,
			  format = "draws_df")) %>%
		select(all_of(param))
	# get HPD bounds
	bounds95 = hdi(p_dat, credMass=0.95)[,1]
	bounds50 = hdi(p_dat, credMass=0.50)[,1]
	bounds = c(-Inf, bounds95[1], bounds50[1], bounds50[2], bounds95[2], Inf)
	true_val = params_dict[param]
	if (param == "prob_MI_baseline"){
		print(param)
		plot_list[[param]] = plot_shadded_hist(p_dat, bounds, approx_bins=50,
			cols=c(unname(args$colors_dict["mi_hpd3"]), 
				unname(args$colors_dict["mi_hpd2"]), 
				unname(args$colors_dict["mi_hpd1"]))) +
		geom_vline(aes(xintercept=!!true_val),lty="11",color='#333333', linewidth=1.5)
	}else{
		plot_list[[param]] = plot_shadded_hist(p_dat, bounds, approx_bins=50) +
		geom_vline(aes(xintercept=!!true_val),lty="11",color='steelblue', linewidth=1.5)
	}
	plot_list[[param]] = plot_list[[param]] +
		xlab(eval(parse(text=plot_params$label[idx])))
	if ((idx - 1)%%args$ncolParams != 0){
		plot_list[[param]] = plot_list[[param]] + 
			scale_y_continuous(breaks=NULL, name=NULL, expand = expansion(mult = c(0, .15)))
	}
}

p = plot_grid(
		plot_grid(p1, 
			plot_grid(NULL, p2, NULL, rel_heights=c(0.15, 0.9, 0.15), ncol=1), rel_widths=c(0.4, 0.75), labels=c('a', 'b')),
			wrap_plots(plot_list, ncol=args$ncolParams) + 
			plot_annotation(tag_levels = list(letters[3:(length(plot_list)+3)], '1')) & 
				theme(
					plot.tag = element_text(face = 'bold', color='#333333'),
					plot.margin = margin(9,9,9,9)),
		ncol=1, rel_heights=c(0.5, (0.3*ceiling(nrow(plot_params)/args$ncolParams))))


ggsave(paste(c('figures/', args$out,'.pdf'), collapse=''), p, width=15, 
	height= 5.45 + (0.3*ceiling(nrow(plot_params)/args$ncolParams) * 10.9))





#p = plot_grid(
#		plot_grid(p1, 
#			plot_grid(NULL, p2,NULL, rel_heights=c(0.15, 0.9, 0.15), ncol=1), rel_widths=c(0.4, 0.75), labels=c('a', 'b')),
#		plot_grid(
#			NULL,
#			wrap_plots(plot_list, ncol=3) + plot_annotation(tag_levels = list(c('c', 'd', 'e', 'f', 'g', 'h'), '1')) & 
#				theme(plot.tag = element_text(face = 'bold', color='#333333')),
#			NULL,
#			ncol=3,
#			rel_widths=c(0.1, 0.9, 0.1)), 
#		ncol=1, rel_heights=c(0.5, 0.6))


#ggsave(paste(c('figures/', args$out,'.pdf'), collapse=''), p, width=15, height=12)
