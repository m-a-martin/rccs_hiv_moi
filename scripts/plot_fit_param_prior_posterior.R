suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(require(cowplot))
suppressMessages(require(patchwork))
suppressMessages(require(HDInterval))
suppressMessages(source('scripts/utils.R'))


p <- arg_parser("plot summary of simulated data")
p <- add_argument(p, "--fit", help="stan fit objects", nargs=1)
p <- add_argument(p, "--plotParams", help="parameters to plot", nargs=1)
p <- add_argument(p, "--nRowPer", help="number of rows per fit item", default=1, nargs=1)
p <- add_argument(p, "--out", help="output figure name", nargs=1)
p <- add_argument(p, "--colors", help="tsv file with color codes", nargs=1, default='config/colors.tsv')
args <- parse_args(p)
#args$colors = 'config/colors.tsv'
#args$plotParams = 'config/tmp.tsv'
#args$out='test'
#args$fit='fit/full_simulation_base_model_.Rds' 
#args$params='simulations/base_simulation_params.tsv' 

cols = read_tsv(args$colors, show_col_types=FALSE)
args$colors_dict = setNames(cols$color, cols$var)

plot_params = read_tsv(args$plotParams, show_col_types=FALSE)

name_conv = c(
	"logit_prob_seq_coeffs[1]" = "logit_seq_prob_lvl_effect")

d = tibble()
ggplot(d) +
  geom_boxplot(aes(
	  x = 0,
      lower = -0.25, 
      upper = 0.25, 
      middle = 0, 
      ymin = -0.5, 
      ymax = 0.5),
    stat = "identity")
ncol = ceiling(nrow(plot_params)/args$nRowPer)

fit = readRDS(args$fit)

# PRIOR DISTRIBUTIONS ARE HARD CODED!
priors = function(x){
	set.seed(1)
	if (x == "logit_prob_seq_baseline"){
		distr=function(x){qnorm(x,0,2)}
		out = distr(c(0.5, 0.025, 0.975, 0.25, 0.75))
	}else if (x == "logit_prob_seq_coeffs"){
		distr=function(x){qnorm(x,0,2)}
		out = distr(c(0.5, 0.025, 0.975, 0.25, 0.75))
	}else if (x == "logit_prob_seq_ind_sd"){
		# not symmetric so we do by simulation
		sim_vals = rcauchy(100000,0,1)
		# only take positive values
		sim_vals = sim_vals[sim_vals>= 0]
		out = c(median(sim_vals), 
			unname(c(hdi(sim_vals, credMass=0.95))),
			unname(c(hdi(sim_vals, credMass=0.50))))
	}else if (x == "logit_prob_mi_baseline"){
		distr=function(x){qnorm(x,0,3.16)}
		# symmetric and unimodal so HPD = percentiles
		out = distr(c(0.5, 0.025, 0.975, 0.25, 0.75))
	}else if (x == "logit_prob_mi_fnr" | x == "logit_prob_mi_fpr" | x == "logit_prob_mi_coeffs"){
		distr=function(x){qnorm(x,0,1)}
		# symmetric and unimodal so HPD = percentiles
		out = distr(c(0.5, 0.025, 0.975, 0.25, 0.75))
	}else if (x == "prob_mi_fnr" | x == "prob_mi_fpr"){
		sim_logit_vals = rnorm(100000,0,1)
		sim_vals = inv_logit(sim_logit_vals)
		out = c(median(sim_vals), 
			unname(c(hdi(sim_vals, credMass=0.95))),
			unname(c(hdi(sim_vals, credMass=0.50))))
	}else if (x == "w"){
		# not symmetric so we do by simulation
		tau = rcauchy(1000000,0,1)
		tau = tau[tau > 0][1:100000]
		xi = rt(1000000,2)
		xi = xi[xi > 0][1:100000]
		sigma = xi * tau
		sim_vals = rnorm(100000, 0, sigma)
		out = c(median(sim_vals), 
			unname(c(hdi(sim_vals, credMass=0.95))),
			unname(c(hdi(sim_vals, credMass=0.50))))
	}
	return(tibble(vals=out, labs=c('median', 'lower95', 'upper95', 'lower50', 'upper50')) %>%
		pivot_wider(names_from='labs', values_from='vals') %>%
		mutate(type="prior"))
}
		

posterior_dat = list()
plots = list()
for (idx in seq(1,nrow(plot_params))){
	param = plot_params$param[idx]
	p_dat = as_tibble(fit$draws(
			variables = c(param),
			  inc_warmup = FALSE,
			  format = "draws_df")) %>%
		select(all_of(param)) %>%
		rename(y := !!param)
	posterior_dat[[param]] = 
		bind_rows(
			tibble(vals = c(
				median(p_dat$y), 
				c(hdi(p_dat$y, credMass=0.95)),
				c(hdi(p_dat$y, credMass=0.50))),
				labs = c('median', 'lower95', 'upper95', 
					'lower50', 'upper50')) %>%
				pivot_wider(names_from=labs, values_from=vals) %>%
				mutate(type="posterior"),
			priors(str_split(param, "\\[", simplify=TRUE)[,1])) %>%
		mutate(
			param = param,
			label = (plot_params %>% filter(param == param))$label[1],
			type=ordered(type, levels=c('posterior', 'prior')))
	plots[[param]] = 
		ggplot(posterior_dat[[param]], aes(x=type, y=median, color=type)) + 
			geom_errorbar( aes(ymin=lower95, ymax=upper95), width=0, linewidth=0.6) +
			geom_errorbar(aes(ymin=lower50, ymax=upper50), width=0, linewidth=2) +
			geom_point(size=3, shape=21, fill='#eaeaea') + 
			scale_color_manual(values=c('#333333', 'peru'), guide='none') +
			ylab(eval(parse(text=plot_params$label[idx]))) +
			coord_flip() + 
			xlab(NULL) +
			gtheme +
			theme(axis.title.x=element_text(size=12))
}

ncol = ceiling(nrow(plot_params)/args$nRowPer)
p = wrap_plots(plots, ncol=ncol)
ggsave(paste(c('figures/', args$out,'.pdf'), collapse=''), p, width=ncol*4.4, 
	height= args$nRowPer*15/4, limitsize = FALSE)


