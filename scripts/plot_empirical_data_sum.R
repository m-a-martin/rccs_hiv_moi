suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(require(tidytree))
suppressMessages(require(ggtree))
suppressMessages(require(ape))
suppressMessages(require(patchwork))
suppressMessages(require(cowplot))
suppressMessages(source('scripts/utils.R'))


plot_rep_tree = function(td_path){
	td_name = str_split(
		tail(t(str_split(td_path, "/", simplify=TRUE)), 1),
		'\\.', simplify=TRUE)[1]
	td_name = paste(c(
		str_split(td_name, "_", simplify=T)[1],
		', ',
		str_split(td_name, "_", simplify=T)[2],
		'-',
		str_split(td_name, "_", simplify=T)[4]), collapse='')
	td = load(td_path)
	p = ggplot(sub_phsc, aes(x,y)) +
		geom_tree(color='#333333') +
		geom_tippoint(color='#333333', aes(fill=subgraph), shape=21, size=3) +
		scale_fill_manual(values=c('#8fbc8f', 'steelblue'), guide="none") +
		xlab('subs/site') +
		scale_y_continuous(breaks=NULL, name=NULL) +
		ggtitle(td_name) +
		gtheme +
		theme(
			axis.line.y = element_blank(),
			axis.line.x = element_line(color='#333333'),
   			panel.grid.major.y = element_blank(),
    		panel.grid.minor.y = element_blank(),
    		panel.border = element_blank(),
    		plot.margin = unit(c(12, 12, 12, 12), "points")) 
	return(p)
}


plot_window_n_d = function(d){
	suppressMessages(require(patchwork))
	window_d = d %>% 
		mutate(window_mid = ( ( window_end + window_start ) / 2 ) ) %>%
		group_by(window_mid) %>%
		summarise(obs = n(), n_obs = sum(n), d_obs = sum(d))
	#mutate(n_obs_p = n_obs/obs, 
	#	d_obs_p = d_obs/n_obs)
	#window_d = 
	#	cbind(window_d, 
	#		as_tibble(
	#				wilson_ci(window_d$n_obs_p, window_d$n_obs, 0.05)) %>%
	#				rename(c('n_obs_ll' = 'll', 'n_obs_ul' = 'ul')),
	#		as_tibble(
	#				wilson_ci(window_d$d_obs_p, window_d$d_obs, 0.05)) %>%
	#				rename(c('d_obs_ll' = 'll', 'd_obs_ul' = 'ul')))
	#p = ggplot(window_d, aes(x=window_mid)) +
	#	geom_line(aes(y=n_obs_p*100), color='#333333', linetype='dashed') +
	#	geom_errorbar(
     #     aes(ymin=n_obs_ll*100,ymax=n_obs_ul*100), color='#333333', width=0) +
	#	geom_line(aes(y=d_obs_p*100), color='maroon4', linetype='dashed') +
	#	geom_errorbar(
    #      aes(ymin=d_obs_ll*100,ymax=d_obs_ul*100), color='#333333', width=0) +
	#	geom_point(data=
	#		window_d %>% select(window_mid, n_obs_p, d_obs_p) %>%
	#			pivot_longer(-window_mid),
	#		aes(x=window_mid, y=value*100, fill=name), color='#333333', size=4, shape=21) +
	#	scale_fill_manual(
	#		values=c('n_obs_p' = '#eaeaea', 'd_obs_p' = 'maroon4'),
	#		labels=c('n_obs_p' = 'sequenced windows', 'd_obs_p' = 'dual infected windows'),
	#		name=NULL) +
	#	xlab('position (nt)') +
	#	ylab('observed (%)') +
	#	ylim(0,100)+
	#	gtheme +
	#	theme(legend.position=c(1,1), legend.justification=c(1,1))
	# todo add spliced proteins!!
	# tat
	# rev
	genome = tibble(
			xmin=c(1, 790, 2085, 5041, 6062, 5559, 5831, 8379, 5970, 8379, 6225, 8797, 9086), 
			xmax=c(634, 2292, 5096, 5619, 6310, 5850, 6045, 8469, 6045, 8653, 8795, 9417, 9719), 
			ymin=c(2.15, 2.15, 0.15, 2.15, 1.15, 0.15, 1.15, 2.15, 0.15, 1.15, 0.15, 2.15, 1.15), 
			ymax=c(2.85, 2.85, 0.85, 2.85, 1.85, 0.85, 1.85, 2.85, 0.85, 1.85, 0.85, 2.85, 1.85), 
			label=c("5' LTR", 'gag', 'pol', 'vif', 'vpu', 'vpr', 'tat',  'tat','rev', 'rev', 'env', 'nef', "3' LTR"))
	genome = genome %>% group_by(label) %>% mutate(n_chunk = n())
	chunked = genome %>% filter(n_chunk == 2)
	connectors = chunked %>% arrange(xmin) %>% group_by(label) %>%
		summarise(
				x1 = xmax[1], 
				x2 = xmax[1]+ 250, 
				x3 = xmin[2] - 250, 
				x4 = xmin[2],
				y1 = ymax[1], 
				y2 = (ymax[1] + ymin[2])/2, 
				y3 = (ymax[1] + ymin[2])/2,
				y4=ymin[2]) %>%
			pivot_longer(-label) %>%
			mutate(idx = substr(name,2,2),
				type=substr(name,1,1)) %>%
			select(label, type, value, idx) %>%
			pivot_wider(names_from=type, values_from=value)
	max_x = max(genome$xmax)
	p0 = ggplot()+
		geom_rect(data=genome, 
			aes(xmin=xmin, ymin=ymin, xmax=xmax, ymax=ymax), color='#333333', fill='#eaeaea') +
		geom_line(data=connectors, 
			aes(x=x, y=y, group=label), color='#707070', linetype='dotted') + 
		geom_text(data = genome %>% filter(n_chunk == 1), 
			aes(x=(xmax + xmin)/2, y=(ymax+ymin)/2, label=label, size=label=='vpu'), color='#333333') +
		geom_text(
			data = 
				connectors %>% filter(idx == 2 | idx == 3) %>%
					group_by(label) %>% summarise(y = y[1], x = sum(x)/n()), 
			aes(x=x, y=y+0.3, label=label), size=4, color='#333333') +
		scale_y_continuous(breaks=NULL, name='ORF') +
		scale_x_continuous(limits=c(0,max_x),expand = c(0.025, 0)) +
		scale_size_manual(values=c('TRUE'=2.75, 'FALSE'=4), guide='none') +
		xlab(NULL) +
		gtheme +
		theme(axis.text.x = element_blank(),
			axis.line.x = element_line(colour = "#333333"),
		    panel.border = element_blank())
	p1 = ggplot(window_d, aes(x=window_mid, y=n_obs)) + 
		geom_bar(stat="identity", fill="#eaeaea", color='#333333') +
		ylab('\nsequenced\nsamples') +
		scale_x_continuous(name=NULL, labels=NULL, limits=c(0,max_x),expand = c(0.025, 0)) +
		xlab(NULL) +
		gtheme +
		theme(axis.text.x = element_blank())
	p2 = ggplot(window_d, aes(x=window_mid, y=d_obs)) + 
		geom_bar(stat="identity", color='#333333', fill=args$colors_dict['multiple_subgraphs']) +
		ylab('\nmultiple subgraph\nsamples') +
		xlab('position (nt)') +
		scale_x_continuous(limits=c(0,max_x),expand = c(0.025, 0)) +
		gtheme 
	#p = plot_grid(p0, p1, p2, ncol=1)
	p = p0 / p1 / p2 + plot_layout(heights=c(0.65, 1, 1)) + 
		plot_annotation(tag_levels = list(c('e', 'f', 'g'), '1')) &
    	theme(plot.tag = element_text(face = 'bold'))
	return(p)
}

plot_cophenetic_dist = function(x, colors_dict){
	x  = x %>% select(id, window_start, id_min_sg_dist) %>%
		filter(!is.na(id_min_sg_dist)) %>%
		unique()
	bounds = c(-Inf, quantile(x$id_min_sg_dist, c(0.025, 0.25, 0.75, 0.975)), Inf)
	p = plot_shadded_hist(x %>% select(id_min_sg_dist), bounds, approx_bins=50,
			cols=c(unname(colors_dict["multiple_subgraphs_1"]), 
				unname(colors_dict["multiple_subgraphs_2"]), 
				unname(colors_dict["multiple_subgraphs_3"]))) +
		xlab('subs/site') +
		geom_vline(xintercept=median(x$id_min_sg_dist), linetype='dashed', color='#333333')
	return(p)
}

p <- arg_parser("plot summary of empirical data")
p <- add_argument(p, "--dat", help="input data file", nargs=1)
p <- add_argument(p, "--out", help="output figure name", nargs=1)
p <- add_argument(p, "--repTree1", help="representative tree data 1", nargs=1)
p <- add_argument(p, "--repTree2", help="representative tree data 2", nargs=1)
p <- add_argument(p, "--colors", help="tsv file with color codes", default='config/colors.tsv', nargs=1)
args <- parse_args(p)
#args$dat = "output/211220_allreads_phsc_all_subgraphs_format.tsv"
#args$repTree1 = 'output/trees/AID2211-fq1_1800_to_2049.Rda'
#args$repTree2 = 'output/trees/AID2642-fq2_1800_to_2049.Rda'
#args$out = 'figure_1'
#args$colors = 'config/colors.tsv'
cols = read_tsv(args$colors, show_col_types=FALSE)
args$colors_dict = setNames(cols$color, cols$var)

dat = read_tsv(args$dat, show_col_types = FALSE) %>% 
	filter(window_type == 'unique')

tab_dat = tabulate_n_d(dat)
sum_dat = summarise_n_d(tab_dat)

p1 = plot_rep_tree(args$repTree1)
p2 = plot_rep_tree(args$repTree2)
p3 = plot_cophenetic_dist(dat, args$colors_dict)


N_win = max(dat$N_windows)
model_pred = tibble(phi = seq(0,1,0.001)) %>%
	mutate(n = N_win*(1-(1-phi)^2), 
		d = N_win*phi^2)

p4 = plot_n_d_scatter(sum_dat, args$colors_dict, grouping="empirical", prediction=model_pred %>% select(n,d)) 
#p4 = plot_n_d(sum_dat)
p5 = plot_window_n_d(tab_dat)


p = plot_grid(
		plot_grid(p1, NULL, p2, rel_widths=c(0.9, 0.05, 0.9), nrow=1, labels=c('a', '', 'b')),
		plot_grid(
			plot_grid(
				NULL,
				p3,
				NULL,
				ncol=1,
				rel_heights=c(0.2,0.9,0.1),
				labels=c('', 'c', '')), 
			NULL, p4, nrow=1, rel_widths=c(0.4,0.1, 0.6), labels=c('', '', 'd')),
		NULL,
		p5,
		ncol=1,
		rel_heights=c(0.3, 0.4, 0.02, 0.4))


ggsave(paste(c('figures/', args$out, '.pdf'), collapse=''), p, width=10, height=14)

