suppressMessages(library(tidyverse))
suppressMessages(library(argparser))
suppressMessages(library(tidytree))
suppressMessages(library(ape))
suppressMessages(source('scripts/utils.R'))


run = function(in_args){
	# this loads alternative args object, 
	# hence use in_args for function args
	suppressMessages(load(in_args$dat))
	window_size = tail_less_head(str_split(names(phyloscanner.trees)[1], '_to_', simplify=TRUE))
	window = 
		paste(c(in_args$windowStart, 
				'_to_',  
				as.character(as.numeric(in_args$windowStart)+window_size)), 
			collapse='')
	phsc = phyloscanner.trees[[window]][['tree']]
	splits = phyloscanner.trees[[window]][['splits.table']] %>%
		mutate(sample = str_split(tip, "_", simplify=TRUE)[,1]) %>%
		filter(sample == in_args$sample) %>%
		select(tip, subgraph) %>%
		rename(label=tip)
	sub_phsc = as.treedata(as_tibble(
		reorder(drop.tip(phsc, setdiff(1:length(phsc$tip.label), 
			which(str_split(phsc$tip.label, '_', simplify=TRUE)[,1] == in_args$sample))))
		) %>%
			left_join(splits, by='label'))
	dir.create(file.path('output', 'trees'), showWarnings = FALSE)
	
	op = file.path('output/trees', paste(c(in_args$sample, '_', window,'.Rda'), collapse=''))
	cat(op)
	save(sub_phsc, file = op)
}


p <- arg_parser("get treedata for a given window and sample")
p <- add_argument(p, "--sample", help="sample of tips to get", nargs=1)
p <- add_argument(p, "--windowStart", help="starting window position", nargs=1)
p <- add_argument(p, "--dat", help="phyloscanner output data file", nargs=1)
args <- parse_args(p)
#args$dat = "data/211220_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd/ptyr298_workspace.rda"
#args$windowStart = "1800"
#args$sample = "AID2642-fq2"
run(args)


