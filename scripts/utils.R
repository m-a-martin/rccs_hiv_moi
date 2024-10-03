

get_hpd = function(x, credMass=0.95){
	suppressMessages(require(HDInterval))
  return(enframe(hdi(x, credMass=credMass)) %>% pivot_wider())
}



read_format_rccs_dat = function(rccs_df){
	rccs_metadata = read_dta(rccs_df) %>%
		mutate(round = as.numeric(str_replace(str_remove(round, "^R0+"), 'S', '.1')))
	# foramt round
	# codes BD copies as -1
	# codes missing copies as NA 
	rccs_metadata = rccs_metadata %>% 
		mutate(
			int_date = 
				as.Date(int_date, 
					format = 
						ifelse(grepl("/", int_date, fixed = TRUE), 
							ifelse(nchar(str_split(int_date, "/", simplify=TRUE)[,3]) == 2, 
								"%d/%m/%y", "%d/%m/%Y"),
							"%Y-%m-%d")),
			copies = gsub(',','',copies),
			copies = if_else(copies == 'BD', -1, as.numeric(copies)),
			log10_copies = case_when(
				is.na(copies) ~ NA,
				copies < 0 ~ copies,
				!is.na(copies) & copies > 0 ~ log10(copies)),
			log10_copies_cat = case_when(
				finalhiv != "P" ~ '-',
				finalhiv == "P" & is.na(log10_copies) ~ 'missing',
				finalhiv == "P" & !is.na(log10_copies) & log10_copies <= log10(200) ~ 'suppressed',
				finalhiv == "P" & !is.na(log10_copies) & log10_copies > log10(200) ~ 
					cut(log10_copies, c(log10(200), 3, 3.5, 4, 4.5, 5, Inf))),
			age_cat = cut(ageyrs, breaks=c(14, 24, 34, 50, 100)),
			comm_type = as_factor(comm_type)) %>%
		# add pre_treatment columns
		left_join(
			rccs_metadata %>% select(study_id, round, arvmed, cuarvmed) %>%
				mutate(arv = if_else(arvmed==1|cuarvmed==1, TRUE, FALSE)) %>%
				filter(!is.na(arv) & arv == TRUE) %>% 
				group_by(study_id) %>%
				summarise(first_arv_round = min(round)),
			by='study_id') %>%
		mutate(pre_treatment = round < first_arv_round | is.na(first_arv_round)) %>%
		filter(!is.na(ageyrs) & ageyrs >= 15 & ageyrs <= 49)
	# bad examples: 
	# H013632  21/06/1995 R001
	# G113713  06/11/2020
	# todo add date formatting!
	return(rccs_metadata)
}



tabulate_n_d = function(x){
	# subgraph specific columns = 
	# sample_subgraph_reads
	# sg
	# dist
	return(x %>%
		group_by(across(-c(id_subgraph_reads, sg, id_min_sg_dist, min_sg_dist))) %>%
		summarise(
			n = 1, 
			n_subgraphs = length(unique(sg)), 
			d = 1*(n_subgraphs > 1), .groups='drop'))
}

summarise_n_d = function(x){
	x = x %>%
		group_by(across(-c(n_subgraphs, n, d, window_start, window_end, window_type))) %>%
		summarise(N_obs = sum(n), 
			MI_obs = sum(d),
			.groups='drop') %>%
		filter(N_obs > 0)
	return(x)
}


tail_less_head = function(x){
	x = as.numeric(x)
	return(tail(x,1)  - head(x,1))
}


inv_logit = function(x){
	return(exp(x)/(1+exp(x)))
}

logit = function(x){
	return(log(x/(1-x)))
}


filter_duplicate_ids = function(x){
	# some IDs are included in multiple phyloscanner runs
	# here we grab just one based on the number of dually infected windows
	# n_cc shouldn't vary much by run, but may be some differences in tree building, etc. 
	x = x %>% right_join(
	x%>% group_by(id, run) %>% 
		summarise(n_dual = sum(n_cc > 1), .groups='drop') %>%
		arrange(-n_dual) %>%
		group_by(id) %>%
		filter(row_number() == 1) %>%
		select(c('id', 'run')),
		by=c("id", "run"))
	return(x)
}



fancy_scientific <- function(l) {
     # turn in to character string in scientific notation
     l <- format(l, scientific = TRUE)
     # quote the part before the exponent to keep all the digits
     l <- gsub("^(.*)e", "'\\1'e", l)
     # turn the 'e+' into plotmath format
     l <- gsub("e", "%*%10^", l)
     # return this as an expression
     parse(text=l)
}


plot_ind_posterior_pred = function(dat, fit, colors_dict, grouping="empirical", guide=FALSE){
	suppressMessages(require(tidyverse))
	fit_dat = as_tibble(fit$draws(
			variables = c("ind_log_prob_mi"),
			  inc_warmup = FALSE,
			  format = "draws_df")) %>%
		mutate(
	  		iter = seq(1,n())) %>%
		pivot_longer(-iter) %>%
		filter(substr(name, 1,15) == 'ind_log_prob_mi')

	# transform
	fit_dat = fit_dat %>% 
		mutate(
			value = exp(value),
			name=gsub('_log', '', name))
		
	# median value
	fit_dat_sum = fit_dat %>%
	  group_by(name) %>% 
	  summarise(median=median(value)) 

	# add bounds, probably cleaner way to do this
	fit_dat_sum = bind_cols(fit_dat_sum, 
	  bind_rows(fit_dat %>% group_by(name) %>%
	  	group_map(~hdi(.x$value, credMass=0.95)))) %>%
	arrange(median) %>%
	mutate(x = seq(1,n()), 
		id = as.numeric(str_replace(str_split(name, "\\[", simplify=TRUE)[,2], '\\]', '')))

	if (grouping == "empirical"){
		fit_dat_sum = fit_dat_sum %>% left_join(
			dat %>% 
				mutate(id = seq(1,n())) %>%
				mutate(group=MI_obs > 0) %>%
				select(id, group),
			by="id")
		labels=c('TRUE'=">= 1 multiple subgraph windows", 'FALSE'="0 multiple subgraph windows")
		highlight_color = unname(colors_dict['multiple_subgraphs'])
	}else if (grouping == "truth"){
		fit_dat_sum = fit_dat_sum %>% left_join(
			dat %>% 
				mutate(group=as.logical(has_MI)) %>%
				select(id, group),
			by="id")
		labels=c('TRUE'="multiply infected", 'FALSE'="singly infected")
		highlight_color = unname(colors_dict['multiply_infected'])
	}
	if (guide == TRUE){
		scale_color = scale_color_manual(values=c("FALSE"="#333333", 
			"TRUE"=highlight_color), 
			labels=labels, name=NULL)
		color_guide = guides(colour = guide_legend(override.aes = list(size=4, alpha = 1), alpha=1, position="inside"))
	}else{
		scale_color = scale_color_manual(values=c("FALSE"="#333333", 
			"TRUE"=highlight_color), 
			labels=labels, guide=NULL)
		color_guide = NULL
	}
	p = ggplot(fit_dat_sum, aes(x=x, y=median, ymin=lower, ymax=upper, color=group)) + 
		geom_errorbar(width=0, color='#adadad') +
		geom_point(shape=21) +
		scale_color + 
		color_guide + 
		scale_x_continuous(breaks=NULL, name="sample") +
		ylab('posterior prob. of\nmultiple infection') +
		gtheme +
		theme(legend.justification.inside=c(0,1))
	return(p)
}

clean_numeric = function(x) ifelse(x < 1 & x > -1, signif(x, 1), round(x))
	
plot_shadded_hist = function(d, bounds, approx_bins=100, cols=c('#eaeaea', '#adadad', "#707070")){
	suppressMessages(require(tidyverse))
	shades = c(cols, rev(cols[1:2]))
	# shades histogram by HPD intervals
	# probably a function that exists for this
	# but fun challenge to rewrite one
	# assumes unimodality
	bw = clean_numeric((max(d[,1]) - min(d[,1]))/approx_bins)
	ll = floor(min(d[,1])/bw)*bw
	ul = ceiling(max(d[,1])/bw)*bw
	breaks = seq(ll, ul, bw)
	d$bin = cut(as_vector(d[,1]), breaks, dig.lab=10)
	sum_d = d %>% 
		group_by(bin) %>%
		summarise(n=n()) %>%
		mutate(
			min = as.numeric(str_replace(str_split(bin, ",", simplify=TRUE)[,1], "\\(", "")),
			max = as.numeric(str_replace(str_split(bin, ",", simplify=TRUE)[,2], "\\]", "")),
			x =  (min + max)/2,
			min_shade_bin = cut(min, bounds, labels=FALSE),
			max_shade_bin = cut(max, bounds, labels=FALSE),
			# below assumes that min_shade_bin <= max_shade_bin
			# five bins :
			# bin 1: 0th to 5th % HPD
			# bin 2: 5th to 25th % HPD
			# bin 3: 25th to 75th % HPD
			# bin 4: 75th to 95th % HPD
			# bin 5: 95th to 100th % HPD
			shade_bin = case_when(
				# if they match, just choose one
				min_shade_bin == max_shade_bin ~ min_shade_bin,
				# if they don't match choose the one that pushes the interval out the most
				(min_shade_bin != max_shade_bin) & 
					min_shade_bin < 3 ~ max_shade_bin,
				(min_shade_bin != max_shade_bin) & 
					min_shade_bin >= 3 ~ min_shade_bin))
	p = ggplot(sum_d, aes(x=x, y=n, fill=as.character(shade_bin))) +
		geom_bar(color='#333333', stat="identity") +
		scale_fill_manual(values=shades,
			guide=NULL) +
		xlab(colnames(d)[1]) +
		scale_y_continuous(breaks=NULL, name='density', expand = expansion(mult = c(0, .15))) +
		gtheme #+
		#theme(panel.grid.major = element_blank(), 
		#	panel.grid.minor=element_blank())
	return(p)
}


eval_labels = function(x){
		out = c()
		for (i in x){
			if (grepl("expression", i) == TRUE){
				out = c(out, eval(parse(text=i)))
			}else{
				out = c(out, i)
			}
		}
		return(out)
}


plot_n_d_scatter = function(d, colors_dict, grouping="empirical", prediction=tibble(x=c(), y=c())){
	suppressMessages(require(patchwork))
	suppressMessages(require(tidyverse))
	set.seed(1)
	d = d %>%
		mutate(N_obs_jitter = N_obs + rnorm(nrow(d), 0,0.2),
			MI_obs_jitter = MI_obs + rnorm(nrow(d), 0, 0.2)) 
	#expression("">=1 multiple subgraph window)
	if (grouping == "empirical"){
		d = d %>% mutate(group = MI_obs > 0)
		labels=c(
			'FALSE'="0 multiple subgraph windows",
			'TRUE'='expression(paste("">=1, " multiple subgraph windows"))')
		highlight_color = unname(colors_dict['multiple_subgraphs'])
	}else if (grouping == "truth"){
		d = d %>% mutate(group = as.logical(has_MI))
		labels=c('FALSE'="singly infected",
			'TRUE'="multiply infected")
		highlight_color = unname(colors_dict['multiply_infected'])
	}
	colnames(prediction) = c('x', 'y')
	p1 = ggplot(d, aes(x=N_obs_jitter, y=MI_obs_jitter, color=group, alpha=group)) +
		geom_point(shape=21)
	if (nrow(prediction) > 0){
		p1 = p1 + geom_line(data=prediction, aes(x=x, y=y), 
			color='#333333', alpha=1, lty="11")
	}
	p1 = p1  +
		scale_color_manual(values=c(
			'TRUE'=highlight_color, 
			'FALSE'='#333333'), labels=eval_labels(labels), name=NULL) +
		guides(colour = guide_legend(override.aes = list(size=4, alpha = 1), alpha=1, position='inside')) +
		scale_alpha_manual(values=c('TRUE'=1, 'FALSE'=0.25), guide='none') + 
		scale_x_continuous(
			#breaks=seq(-1, max(d$N_obs)+1, 1), 
			limits=c(-1, max(d$N_obs)+1),
			name=expression(paste("sequenced windows (", N[i]^{obs}, ")"))) +
		scale_y_continuous(
			#breaks=seq(-1, max(d$N_obs)+1, 1), 
			limits=c(-1, max(d$N_obs)+1),
			name=expression(atop("multiple subgraph", paste("windows (MI" [i]^{obs}, ")")))) +
		gtheme +
		#theme(aspect.ratio=1) +
	 	theme(plot.margin = margin(0, 0, 0, 0, "pt"),
	 		legend.justification.inside=c(0,1))

	p2 = ggplot(d %>% group_by(MI_obs) %>% summarise(n=n()), aes(x=MI_obs, y=n)) +
		geom_bar(stat="identity", color='#333333', fill='#eaeaea') +
		coord_flip() + 
		scale_x_continuous(breaks=NULL, name=NULL, limits=c(-1, 30)) +
		scale_y_continuous(breaks=NULL, name=expression(log[10](density)), trans = 'log10',
			expand = expansion(mult = c(0, .1))) +
		gtheme +
	  theme(
	  	plot.margin = margin(0, 0, 0, 0, "pt"),
	    panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	    panel.background = element_blank())

	p3 = ggplot(d %>% group_by(N_obs) %>% summarise(n=n()), aes(x=N_obs, y=n)) +
		geom_bar(stat="identity", color='#333333', fill='#eaeaea') +
		scale_x_continuous(breaks=NULL, name=NULL, limits=c(-1, 30)) +
		scale_y_continuous(breaks=NULL, name='density',
			expand = expansion(mult = c(0, .1))) +
		gtheme +
	  theme(plot.margin = margin(0, 0, 0, 0, "pt"),
	    panel.grid.major = element_blank(),
	    panel.grid.minor = element_blank(),
	    panel.border = element_blank(),
	    panel.background = element_blank())

	p = (p3 + 
			( plot_spacer() + 
				theme(plot.margin = margin(0, 0, 0, 0, "pt")) ) + 
			plot_layout(widths=c(1,0.5))) / 
		(p1 + p2 + plot_layout(widths=c(1,0.5))) + plot_layout(heights=c(0.5,1))

	return(p)
}


gtheme = theme_bw() + 
	theme(text=element_text(color='#333333'),
			plot.title = element_text(hjust = 0.5, size=14),
			axis.text=element_text(color='#707070', size=12),
			legend.text = element_text(color='#333333', size=12),
			legend.title = element_text(color='#333333', size=12),
			axis.title=element_text(size=16),
			axis.ticks = element_line(color = "#707070"),
			panel.grid.major = element_line('#eaeaea', 0.5),
		  panel.grid.minor = element_line('#eaeaea', 0.5),
		  legend.background = element_rect(fill='transparent'),
		  strip.background = element_blank(), 
		  strip.text = element_text(size=14))
