suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(require(cowplot))
suppressMessages(require(bayesplot))
suppressMessages(require(posterior, include.only=c("subset_draws")))
suppressMessages(source('scripts/utils.R'))


p <- arg_parser("plot summary of simulated data")
p <- add_argument(p, "--fit", help="stan fit object", nargs=1)
p <- add_argument(p, "--plotParams", help="parameters to plot", nargs=1)
p <- add_argument(p, "--samples", help="number of posterior samples per chain to plot", 
	default=250, nargs=1, type="integer")
p <- add_argument(p, "--out", help="output figure name", nargs=1)
args <- parse_args(p)
#args$fit = 'fit/base_simulation_base_model.Rds'
#args$plotParams = 'config/base_simulation_plot_params.tsv'
plot_params = read_tsv(args$plotParam, show_col_types=FALSE)

fit = readRDS(args$fit)
po <- fit$draws(
  variables = c("lp__",plot_params$param),
  inc_warmup = FALSE,
  format = "draws_df")

# sample iterations
po = posterior::subset_draws(po, iteration=sample.int(max(po$.iteration), args$samples, replace=FALSE))

colnames(po) = c('log(posterior)', plot_params$simple_label, colnames(po)[(ncol(po)-2):ncol(po)])

# relabel column names
bayesplot_theme_set(
	theme_minimal() +
		theme( 
			panel.grid.major = element_line('#eaeaea', 0.5),
		 	panel.grid.minor = element_line('#eaeaea', 0.5),
		 	text=element_text(color='#333333'),
		 	axis.text=element_text(color='#707070', size=8)))

color_scheme_set(c('#eaeaea', '#eaeaea', '#eaeaea', '#333333', '#333333', '#333333'))


p <- bayesplot::mcmc_pairs(po, 
    off_diag_args = list(size = 0.3, alpha = 0.3))

ggsave(paste(c('figures/', args$out,'.pdf'), collapse=''), p, 
	width=nrow(plot_params)*2.6, height=nrow(plot_params)*2.6,
	limitsize = FALSE)


