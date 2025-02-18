library(tidyverse)
suppressMessages(source('scripts/utils.R'))


dat = read_tsv('output/211220_allreads_phsc_all_subgraphs_dists.tsv', show_col_types=FALSE)
col = (read_tsv('config/colors.tsv', show_col_types=FALSE) %>% filter(var == "multiple_subgraphs"))$color

bins = c(-Inf, c(0, 1, 0.05), Inf)

sum_dat = dat %>%
	mutate(
		bin_end = as.numeric(gsub("]", "", str_split(bin, ",", simplify=TRUE)[,2])),
		coarse_bin_end = round(bin_end * 100)/100) %>%
	group_by(coarse_bin_end, type) %>%
	summarise(n=sum(n)) %>%
	group_by(type) %>%
	mutate(p = n / sum(n),
		type = ordered(type, levels=c('within_sg', 'between_sg'))) 

p = ggplot(sum_dat, aes(x=coarse_bin_end, y=p, fill=type, group=type)) +
	geom_bar(position="identity", stat="identity", alpha=0.75, color='#333333') +
	xlab('subs/site') +
	ylab('density') +
	scale_y_continuous(breaks=NULL, expand = expansion(mult = c(0, .15))) +
	scale_fill_manual(values=c('within_sg' = '#eaeaea', 'between_sg' = col),
		labels=c('within_sg' = 'within subgraph', 'between_sg' = "between subgraph"),
		name=NULL) +
	guides(fill = guide_legend(position = "inside")) + 
	gtheme +
	theme(legend.position.inside = c(1, 1),
		legend.justification.inside=c(1,1))


ggsave('figures/empirical_cophenetic_dists.pdf', width=8.4, height=4.8)
