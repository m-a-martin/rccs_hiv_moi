suppressMessages(library(tidyverse))
suppressMessages(library(ape))
#suppressMessages(library(igraph))
suppressMessages(library(phytools))
# needed for distRoot function
# todo reduce dependencies
suppressMessages(library(adephylo))
suppressMessages(source('scripts/utils.R'))


getMRCA_robust = function(tree, tips){
	if (length(tips) == 1){
		return(tips[1])
	}else{
		return(getMRCA(tree, tips))
	}
}


get_mrca_bl = function(wdat, tips){
	if (length(tips) == 1){
		return(0)
	}else{
		# get just these tips from tree 
		subtree = keep.tip(wdat$tree, tips)
		mrca_branches = subtree$edge[which(subtree$edge[,1] == getMRCA_robust(subtree, tips)),]
		mrca_branches = cbind(mrca_branches, rep(NA, nrow(mrca_branches)))
		for (i in 1:nrow(mrca_branches)){
			mrca_branches[i,3] = nodeheight(subtree, mrca_branches[i,2])
		}
		return(sum(sort(mrca_branches[,3])[1:2]))
	}
}


#### MAIN CODE STARTS HERE ####
in_args = commandArgs(trailingOnly=TRUE)

print(in_args)

#in_args = c(
#	"test.tsv",
#	"data/211220_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd/ptyr298_workspace.rda",
#	"data/211220_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd/ptyr2_workspace.rda")


all_mrca_dists = tibble(
  sample = character(),
  dist = numeric(),
  run = character(),
  window_start = numeric(),
  window_end = numeric())


for (i in 2:length(in_args)){
	print(i)
	print(in_args[i])
	load(file=in_args[i])
	run = str_split(rev(str_split(in_args[i], '/', simplify=T))[1], "_", simplify=T)[,1]
	# iterate over .rda files, calc pi for each, add window label, combine across rda files
	for (window in sort(names(phyloscanner.trees))){
		window_start = as.numeric(str_split(window, "_", simplify=T)[,1])
		window_end = as.numeric(str_split(window, '_', simplify=T)[,3])
		# get a list of all ids in the tree for this window
		window_ids = phyloscanner.trees[[window]]$bl.report %>%
			filter(kept & !grepl("^CNTRL", tip) & !grepl("REF", tip))
		if (nrow(window_ids) > 0){
			window_ids = window_ids %>%
				mutate(sample = str_split(tip, "_", simplify=TRUE)[,1])
					mrca_dists = bind_rows(window_ids %>% group_by(sample) %>%
			group_map(~tibble(sample=.y[[1]], dist=get_mrca_bl(phyloscanner.trees[[window]], .x$tip)))) %>%
			mutate(run=run, window_start = window_start, window_end = window_end)
			all_mrca_dists = bind_rows(all_mrca_dists, mrca_dists)
		}
	}
}

write_tsv(all_mrca_dists, in_args[1])


