suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(require(cowplot))
suppressMessages(require(patchwork))
suppressMessages(require(HDInterval))
suppressMessages(source('scripts/utils.R'))


p <- arg_parser("plot summary of simulated data")
p <- add_argument(p, "--fit", help="stan fit objects", nargs=1)
p <- add_argument(p, "--params", help="simulation parameters", nargs=1)
p <- add_argument(p, "--plotParams", help="parameters to plot", nargs=1)
p <- add_argument(p, "--nRowPer", help="number of rows per fit item", default=1, nargs=1)
p <- add_argument(p, "--out", help="output figure name", nargs=1)
p <- add_argument(p, "--colors", help="tsv file with color codes", nargs=1, default='config/colors.tsv')
args <- parse_args(p)
#args$colors = 'config/colors.tsv'
#args$plotParams = 'config/tmp.tsv'
#args$out='test'
#args$fit='fit/full_simulation_base_model.Rds' 
#args$params='simulations/base_simulation_params.tsv' 

cols = read_tsv(args$colors, show_col_types=FALSE)
args$colors_dict = setNames(cols$color, cols$var)

plot_params = read_tsv(args$plotParams, show_col_types=FALSE)

name_conv = c(
	"logit_prob_seq_coeffs[1]" = "logit_seq_prob_lvl_effect")

#args$fit = 'fit/sensitivity/**.Rda'
#args$params = 'simulations/sensitivity/*delta*_params.tsv'
fit_paths = sort(Sys.glob(args$fit))

if (!is.na(args$param)){
	params_paths = sort(Sys.glob(args$params))
}

ncol = ceiling(nrow(plot_params)/args$nRowPer)

plot_list = list()
for (jdx in 1:length(fit_paths)){
	fit_path = fit_paths[jdx]
	fit_name = str_split(fit_path, "/", simplify=TRUE)
	fit_name = str_replace(fit_name[,ncol(fit_name)], ".Rda", "")
	fit = readRDS(fit_path)
	if (!is.na(args$param)){
		params = read_tsv(params_paths[jdx], show_col_types=FALSE) %>% filter(arg != 'out')
		params_dict = setNames(as.numeric(params$value), params$arg)
	}
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
		# hacky way to clean the bounds
		# todo revisit this
		if (bounds95[1] == bounds50[1]){
			bounds[3] = bounds[3]*1.00001
		}
		if (bounds95[2] == bounds50[2]){
			bounds[4] = bounds[4]*0.99999
		}
		if (!is.na(args$param)){
			if (param %in% names(name_conv)){
				true_val = params_dict[[name_conv[param]]]
			}else{
				true_val = params_dict[[param]]
			}
		}
		if (param == "prob_MI_baseline"){
			plot_list[[length(plot_list)+1]] = plot_shadded_hist(p_dat, bounds, approx_bins=50,
				cols=c(unname(args$colors_dict["mi_hpd3"]), 
					unname(args$colors_dict["mi_hpd2"]), 
					unname(args$colors_dict["mi_hpd1"])))
			if (!is.na(args$param)){
				plot_list[[length(plot_list)]] = plot_list[[length(plot_list)]] + 
					geom_vline(aes(xintercept=!!true_val),lty="11",color='#333333', linewidth=1.5)
			}
		}else{
			plot_list[[length(plot_list)+1]] = plot_shadded_hist(p_dat, bounds, approx_bins=50)
			if (!is.na(args$param)){
				plot_list[[length(plot_list)]] = plot_list[[length(plot_list)]] + 
					geom_vline(aes(xintercept=!!true_val),lty="11",color='steelblue', linewidth=1.5)
			}
		}
		plot_list[[length(plot_list)]] = plot_list[[length(plot_list)]] +
			xlab(eval(parse(text=plot_params$label[idx])))
		if (idx%%ncol != 1){
			plot_list[[length(plot_list)]] = plot_list[[length(plot_list)]] + 
				scale_y_continuous(breaks=NULL, name='', expand = expansion(mult = c(0, .15)))
		}
		if (jdx != length(fit_paths)){
			plot_list[[length(plot_list)]] = plot_list[[length(plot_list)]] + 
				xlab('')
		}
	}
}


p = wrap_plots(plot_list, ncol=ncol)
ggsave(paste(c('figures/', args$out,'.pdf'), collapse=''), p, width=nrow(plot_params)*20/6, 
	height= length(fit_paths)*args$nRowPer*15/4, limitsize = FALSE)

