suppressMessages(library(tidyverse))
suppressMessages(library(ape))
#suppressMessages(library(igraph))
suppressMessages(library(phytools))
# needed for distRoot function
# todo reduce dependencies
suppressMessages(library(adephylo))
suppressMessages(source('scripts/utils.R'))


get_cophenetic_dists = function(x, keep){
    if (length(keep) > 1){
        x = keep.tip(x, keep)
        c = cophenetic(x)
        c[upper.tri(c, diag = FALSE)] = NA
        return(c)
    }else{return(NA)}
}


in_args = commandArgs(trailingOnly=TRUE)

print(in_args)

#in_args = c(
#   "",
#    'output/211220_allreads_phsc_all_subgraphs_format_par.tsv',
#    "data/211220_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd/ptyr100_workspace.rda",
 #   "data/211220_phsc_phscrelationships_02_05_30_min_read_null_max_read_posthoccount_im_mrca_fixpd/ptyr2_workspace.rda")

metadata = read_tsv(in_args[2], show_col_types=FALSE) %>%
    filter(window_type == "unique") %>%
    select(id, run, window_start, window_end) %>%
    unique()

subgraph_dists = tibble(
    type = character(),
    run = character(),
    window_start = numeric(),
    window_end = numeric(),
    bin = character(),
    n = numeric())

bins = c(-Inf, seq(0, 1, 0.001), Inf)

for (i in 3:length(in_args)){
    print(in_args[i])
    load(file=in_args[i])
    i_run = str_split(rev(str_split(in_args[i], '/', simplify=T))[1], "_", simplify=T)[,1]
    run_metadata = metadata %>% filter(run == i_run)
    run_windows = run_metadata %>% select(window_start, window_end) %>% unique()
    keep_id = run_metadata %>% select(id) %>% unique()
    if (nrow(run_windows) > 0){
        for (idx in 1:nrow(run_windows)){
            window_start = as.numeric(run_windows$window_start[idx])
            window_end = as.numeric(run_windows$window_end[idx])
            window = paste(as.character(window_start), 
                as.character(window_end), sep='_to_')
            # get a list of all ids in the tree for this window
            window_subgraphs = phyloscanner.trees[[window]]$duals.info %>% 
                left_join(phyloscanner.trees[[window]]$bl.report, by=c('tip.name'='tip')) %>%
                filter(!grepl("^CNTRL", tip.name) & kept == TRUE) %>%
                select(tip.name, split.ids)         
            if (nrow(window_subgraphs) > 0){
                window_subgraphs = window_subgraphs %>%
                    mutate(
                        id = str_split(tip.name, "_", simplify=T)[,1]) %>%
                    inner_join(keep_id, by='id')
                if (nrow(window_subgraphs) > 1){
                    # for each id get cophenetic distances labelled by w/in or b/w subgraph
                    window_subgraph_dists = bind_rows(window_subgraphs %>%
                        group_by(id) %>%
                        group_map(~
                            as_tibble(get_cophenetic_dists(phyloscanner.trees[[window]]$tree, .x$tip.name))  %>%
                                mutate(tip.name1 = colnames(.)) %>%
                                pivot_longer(-tip.name1, names_to='tip.name2', values_to='dist'))) %>%
                        filter(!is.na(dist) & tip.name1 != tip.name2)
                    if (nrow(window_subgraph_dists) > 0){
                        subgraph_dists = bind_rows(
                            subgraph_dists,
                            window_subgraph_dists %>%
                                left_join(
                                    window_subgraphs %>%
                                        rename(tip.name1 = tip.name, split.ids1 = split.ids),
                                    by=c('tip.name1')) %>%
                                left_join(
                                    window_subgraphs %>%
                                        rename(tip.name2 = tip.name, split.ids2 = split.ids),
                                    by=c('id', 'tip.name2')) %>%
                                mutate(
                                    type = if_else(split.ids1 == split.ids2, 'within_sg', 'between_sg'),
                                    bin=cut(dist, breaks=bins)) %>%
                                group_by(type, bin) %>%
                                summarise(
                                    n=n(), .groups='drop') %>%
                                mutate(window_start = window_start, window_end = window_end, run=i_run))
                    }
                }
            }
        }
    }
}


subgraph_dists = subgraph_dists %>%
    group_by(window_start, window_end, bin, type) %>%
    summarise(n = sum(n), .groups='drop') 


write_tsv(subgraph_dists, in_args[1])

