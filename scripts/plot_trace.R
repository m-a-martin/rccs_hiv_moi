suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(require(cowplot))
suppressMessages(source('scripts/utils.R'))


p <- arg_parser("plot summary of simulated data")
p <- add_argument(p, "--fit", help="stan fit object", nargs=1)
p <- add_argument(p, "--plotParams", help="parameters to plot", nargs=1)
p <- add_argument(p, "--out", help="output figure name", nargs=1)
p <- add_argument(p, "--nRow", help="number of rows", default=1, nargs=1)
args <- parse_args(p)
#args$fit = 'fit/base_simulation_base_model.Rda'
#args$plotParams = 'config/base_simulation_plot_params.tsv'
plot_params = read_tsv(args$plotParam, show_col_types=FALSE)
fit = readRDS(args$fit)

po <- fit$draws(
  variables = c("lp__",plot_params$param),
  inc_warmup = FALSE,
  format = "draws_df")

po2 = po %>% mutate(chain = as.character(po$.chain), iter = po$.iteration) %>%
	pivot_longer(-c(chain, iter)) %>%
	filter(name %in% c("lp__",plot_params$param)) %>%
	mutate(name = 
		ordered(
			name,
			levels=c('lp__', plot_params$param),
			labels=c('log(posterior)', eval_labels(plot_params$label))))

p = ggplot(po2, aes(x=iter, y=value, group=chain, color=chain)) +
	geom_line(alpha=0.5, linewidth=0.5) + 
	facet_wrap(~name, nrow=args$nRow, scales='free_y',labeller = label_parsed) +
	#facet_grid(rows=vars(name), scales='free_y',labeller = label_parsed) +
	xlab('iteration') +
	ylab('value') +
	#scale_color_manual(values=c("#1c3448", "#315b7d", "#4682b4", "#7da7ca")) +
	scale_color_manual(values=c('#232323', '#333333', '#707070', '#adadad')) +
	gtheme


ggsave(paste(c('figures/', args$out,'.pdf'), collapse=''), p, 
	width=18, height=4.8*args$nRow, limitsize = FALSE)

