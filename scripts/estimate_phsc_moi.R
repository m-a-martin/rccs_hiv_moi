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


getHeight_robust = function(tree, obj){
	if (suppressWarnings(is.na(as.numeric(obj)))){
		return(suppressWarnings(distRoot(tree, obj)))
	}else{
		return(nodeheight(tree,as.numeric(obj)))
	}
}


get_sg_dist = function(wdat, id, split_ids){
	mrcas = wdat$duals.info %>% 
		filter(split.ids %in% split_ids) %>%
		group_by(split.ids) %>%
		mutate(
			mrca = as.character(getMRCA_robust(wdat$tree, tip.name)),
			height = getHeight_robust(wdat$tree, mrca[1])) %>%
		ungroup() %>%
		mutate(
			all_mrca = getMRCA_robust(wdat$tree, tip.name),
			all_mrca_height = getHeight_robust(wdat$tree, all_mrca[1])) %>%
		select(split.ids, mrca, height, all_mrca, all_mrca_height) %>%
		unique()
	out = tibble(id=id, dist=sum(mrcas$height - mrcas$all_mrca_height))
	return(out)
}


in_args = commandArgs(trailingOnly=TRUE)

print(in_args)

#in_args = c(
#	"",
#	"data/211220_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd/ptyr298_workspace.rda",
#	"data/211220_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd/ptyr2_workspace.rda")


all_subgraphs = tibble(
  split.ids = character(),
  id = character(),
  id_subgraph_reads = numeric(),
	window_start = numeric(),
	window_end = numeric(),
	run = character(),
	id_min_sg_dist = numeric(),
	min_sg_dist = numeric())

ids = tibble(
	run = character(),
	id = character())

refs = tibble(
	run = character(),
	ref = character())

for (i in 2:length(in_args)){
	print(in_args[i])
	load(file=in_args[i])
	run = str_split(rev(str_split(in_args[i], '/', simplify=T))[1], "_", simplify=T)[,1]
	# iterate over .rda files, calc pi for each, add window label, combine across rda files
	for (window in sort(names(phyloscanner.trees))){
		window_start = as.numeric(str_split(window, "_", simplify=T)[,1])
		window_end = as.numeric(str_split(window, '_', simplify=T)[,3])
		# get a list of all ids in the tree for this window
		window_tips = phyloscanner.trees[[window]]$original.tip.labels
		window_ids = window_tips[!grepl("^REF", window_tips) & !(grepl("^CNTRL", window_tips))]
		if (length(window_ids) > 0){
			window_ids = unique(str_split(window_ids, "_",
				simplify=TRUE)[,1])
			ids = bind_rows(
				ids,
				tibble(id=window_ids) %>%
					mutate(run = run)) %>%
				unique()
		}
		window_refs = window_tips[grepl("^REF", window_tips)]
		if (length(window_refs) > 0){
			refs = bind_rows(
				refs,
				tibble(ref=window_refs) %>%
					mutate(run = run)) %>%
				unique()
		}
		window_subgraphs = phyloscanner.trees[[window]]$duals.info %>% 
			left_join(phyloscanner.trees[[window]]$bl.report, by=c('tip.name'='tip')) %>%
			filter(kept & !grepl("^CNTRL", tip.name)) 
		if (nrow(window_subgraphs) > 0){
			window_subgraphs = window_subgraphs %>%
				mutate(
					id = str_split(tip.name, "_", simplify=T)[,1],
					tip_reads = as.numeric(str_split(tip.name, "_", simplify=T)[,5])) %>%
				# duals.info tibble does not account for contamination filter
				# so need to remerge with tips to get those tips that have been labelled as contaminants
				group_by(split.ids, id) %>% 
				summarise(id_subgraph_reads = sum(tip_reads), .groups='drop')
			# get cophenetic distance between MRCA of two largeset subgraphs in each sample
			duals = window_subgraphs %>% 
				group_by(id) %>% filter(n() > 1) %>%
				ungroup() %>%
				arrange(-id_subgraph_reads) %>% 
				group_by(id) %>% 
				slice_head(n=2)
			if (nrow(duals) > 0){
				dists = bind_rows(duals %>% group_by(id) %>%
					group_map(~get_sg_dist(phyloscanner.trees[[window]], .y[[1]], .x$split.ids))) %>%
					rename("id_min_sg_dist" = "dist")
				window_subgraphs = window_subgraphs %>% left_join(dists, by='id')
			}
			# for each subgraph, need to get minimum distance to another subgraph in the tree
			# phyloscanner does this, but does it for each "subgraph" estimated in the non-pruned tree
			# need to map these subgraphs back to the split.ids estimated in the pruned tree
			# we do this via tip assignments
			split_subgraph_map = phyloscanner.trees[[window]]$splits.table %>% select(tip, subgraph) %>%
					inner_join(
						phyloscanner.trees[[window]]$duals.info %>% select(tip.name, split.ids),
						by=c('tip'='tip.name'))	%>%
					select(subgraph, split.ids) %>%
					unique()
			# next get minimum distances for each subgraph to all other subgraphs in the tree
			# treat each subgraph as both a descendant and a parent
			# then merge withi splits and get minimum for each split
			min_split_dists = bind_rows(
				phyloscanner.trees[[window]]$classification$collapsed %>% 
					select(unique.split, length) %>%
					rename("subgraph"="unique.split"),
				phyloscanner.trees[[window]]$classification$collapsed %>% 
					select(parent.split, length) %>%
					rename("subgraph"="parent.split")) %>%
				left_join(split_subgraph_map, by="subgraph") %>%
				group_by(split.ids) %>%
				summarise(min_sg_dist=min(length), .groups="keep") %>%
				slice(1) %>%
				ungroup() 
			window_subgraphs = window_subgraphs %>% left_join(min_split_dists, by='split.ids')
			all_subgraphs = bind_rows(
				all_subgraphs,
				window_subgraphs %>%
					mutate(
						window_start = window_start, 
						window_end = window_end,
						run = run))
		}
	}
}


all_subgraphs = all_subgraphs %>%
	arrange(-id_subgraph_reads) 
	
write_tsv(all_subgraphs, in_args[1])

write_tsv(ids, str_replace(in_args[1], '.tsv', '_ids.tsv'))
write_tsv(refs, str_replace(in_args[1], '.tsv', '_refs.tsv'))

