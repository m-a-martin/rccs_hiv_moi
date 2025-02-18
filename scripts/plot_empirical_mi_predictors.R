suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(require(cowplot))
suppressMessages(require(patchwork))
suppressMessages(require(HDInterval))
suppressMessages(source('scripts/utils.R'))


get_hpd = function(x, credMass=0.95){
  return(enframe(hdi(x, credMass=credMass)) %>% pivot_wider())
}


calc_hpd_hist_bounds = function(d){
    bounds95 = hdi(d, credMass=0.95)
    bounds50 = hdi(d, credMass=0.50)
    return(c(-Inf, bounds95[1], bounds50[1], bounds50[2], bounds95[2], Inf))
}


summarise_draws = function(fit_draws){
    fit_sum = 
        bind_rows(
            fit_draws %>% 
                group_by(name) %>% 
                group_map(~get_hpd(.x$value) %>% 
                  mutate(name=.y$name))) %>%
        left_join(
            bind_rows(
                fit_draws %>% 
                    group_by(name) %>% 
                    group_map(~get_hpd(.x$value, credMass=0.5) %>% 
                      mutate(name=.y$name))) %>%
                rename(c('lower50' = 'lower', 'upper50' = 'upper')),
            by='name') %>%
        left_join(
            fit_draws %>% group_by(name) %>%
                    summarise(
                  median = median(value)),
              by='name')
    return(fit_sum)
}

# check out in-migrant coefficients!



plot_coeffs = function(fit_draws, mi_dm, colors_dict){
    # subset draws to just the columns we want  
    fit_draws = fit_draws[,grepl('logit_prob_mi_coeffs', colnames(fit_draws))] %>%
        rename_with(~gsub('logit_prob_mi_coeffs\\[', '', gsub('\\]', '', .))) %>%
        mutate(cit = seq(1,n())) %>%
        pivot_longer(-cit) %>%
        mutate(
            name = mi_dm$col[as.numeric(name)],
            col1 = str_split(name, '\\*', simplify=TRUE)[,1],
            col2 = str_split(name, '\\*', simplify=TRUE)[,2])
    # need to group interaction terms to simplify figure
    # assume all int_terms share grouping variable (here, community type)
    int_terms = unique(str_split(mi_dm$col[grepl('\\*', mi_dm$col)], '\\*', simplify=TRUE)[,1])
    # define output
    int_coeffs = tibble()
    comms = unique(fit_draws$col2)
    comms = comms[comms != ""]
    for (comm in comms){
        for (var in int_terms){
            int_coeffs = bind_rows(
                int_coeffs, 
                fit_draws %>% filter( name == var | name == paste(var, comm, sep='*')) %>%
                    group_by(cit) %>%
                    summarise(value = sum(value), name=paste(var, comm, sep='*'), .groups='drop'))
        }
    }
    # drop original interaction coefficients and row bind tibbles
    fit_draws = bind_rows(
        fit_draws %>%
            filter(!(col1 %in% int_terms)),
        int_coeffs)
     # summarise draws
    fit_sum = summarise_draws(fit_draws %>% select(-col1, -col2))
    # get sorted y values for each value
    y_vals = mi_dm %>%
        mutate(cat = ordered(
                case_when(
                    grepl('comm_typeinland', col) ~ 'inland',
                    grepl('comm_typefishing', col) ~ 'fishing',
                    TRUE ~ ''),
                levels=c("", 'inland', 'fishing'))) %>%
        arrange(cat, 
            !grepl('sequencing_technology', col), 
            !grepl('round', col), 
            !(col == 'comm_typeinland' | col == 'comm_typefishing'),
            col) %>% 
        filter(col %in% unique(fit_sum$name)) %>% 
        mutate(y = seq(1,n())) %>%
        select(col, y)
    # categorize coefficients and clean up labels,
    # I apologize for verbosity here 
    fit_sum = fit_sum %>%
        left_join(y_vals,
            by=c('name' = 'col')) %>%
        mutate(cat = ordered(
                case_when(
                    grepl('comm_typeinland', name) ~ 'inland',
                    grepl('comm_typefishing', name) ~ 'fishing',
                    TRUE ~ ''),
                levels=c("", 'inland', 'fishing'))) %>%
        mutate(col = str_split(name, '\\*', simplify=TRUE)[,1]) %>%
        mutate(
            col = gsub('round14', '2010 survey', col),
            col = gsub('round15', '2012 survey', col),
            col = gsub('round16', '2014 survey', col),
            col = gsub('round17', '2015 survey', col),
            col = gsub('round18', '2017 survey', col),
            col = gsub('round19', '2019 survey', col),
            col = gsub('barworkerTRUE', 'sex & bar/rest. work', col),
            col = gsub('barworkerFALSE', 'no sex & bar/rest. work', col),
            col = gsub('in_migrantFALSE', 'non-migrant', col),
            col = gsub('male_circumcisionTRUE', 'circumcised men', col),
            col = gsub('male_circumcisionFALSE', 'uncircumcised men', col),
            col = gsub("TRUE", "", col),
            col = gsub('sexM', 'men ', col),
            col = gsub('sexF', 'women ', col),
            col = gsub('plhiv_sexpever_std', 'lifetime sex partners', col),
            col = gsub('comm_typeinland', 'baseline', col),
            col = gsub('comm_typefishing', 'baseline', col),
            col = gsub('fishing', '', col),
            col = gsub('inland', '', col),
            col = gsub('sequencing_technology', '', col),
            col = gsub('age_cat_coarse', 'age ', col),
            col = gsub('_', ' ', col),
            col = gsub(':', '', col),
            col=if_else(col == "", 'intercept', col)) %>%
        #group_by(cat) %>%
        arrange(-y) %>%
        mutate(col = ordered(col, levels=unique(col)))
    rel_heights = (fit_sum %>% group_by(cat) %>% summarise(n=n()))$n
    xspan = max(fit_sum$upper) - min(fit_sum$lower)
    xlim = c(min(fit_sum$lower)-xspan*0.025, max(fit_sum$upper)+xspan*0.025)
    plotlist = list()
    for (i in rev(unique(fit_sum$cat))){
        if (i == ""){
            color = '#333333'
        }else{color = colors_dict[i]}
        plotlist[[length(plotlist) + 1]] = ggplot(fit_sum %>% filter(cat == i), 
                aes(x=median, y=col)) +
            geom_vline(aes(xintercept=0), color='#adadad',linewidth=0.5) +
            geom_errorbar(aes(xmin=lower, xmax=upper), color=color, width=0, linewidth=0.6) +
            geom_errorbar(aes(xmin=lower50, xmax=upper50), color=color, width=0, linewidth=2) +
            geom_point(fill='#eaeaea', color=color, shape=21, size=3) +
            scale_y_discrete(name=i) +
            scale_x_continuous(name=NULL, limits=xlim) +
            gtheme +
            theme(
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                panel.grid.minor.y = element_blank(),
                panel.grid.major.y = element_blank(),
                panel.border = element_blank(),
                axis.line.y = element_line(colour = "#333333"),
                plot.margin = margin(0, 0, 0, 0, "pt"))
    }
    plotlist[[length(plotlist)]] = plotlist[[length(plotlist)]] + 
        scale_x_continuous(name='coefficient', limits=xlim) + 
        theme(
            axis.line.x = element_line(colour = "#333333"),
            axis.text.x=element_text(color='#707070', size=12),
            axis.ticks.x = element_line(color = "#707070"))
    return(wrap_plots(plotlist, ncol=1) + 
                theme(plot.margin = margin(0, 0, 0, 0, "pt")) +
                plot_layout(heights=rel_heights))
}


plot_comm_type = function(draws, dat, metadata, colors_dict, strata_cols = c('sex', 'age_cat_coarse', 'comm_type'), approx_bins=50){
    # implicitly assumes no rows were dropped when building design matrix
    # fine for sex, age_cat_coarse, and comm_type, may be not for other design matrices
    uniq_dat = dat %>% select(all_of(strata_cols)) %>%
                mutate(idx = as.character(seq(1,n()))) %>%
                group_by_at(strata_cols) %>%
                slice(1)
    prop_mi = draws[grepl('prob_mi\\[', colnames(draws)) & !grepl('logit', colnames(draws))] %>%
        rename_with(~gsub("prob_mi\\[", "", gsub("\\]", "", .))) %>%
        mutate(cit = seq(1,n())) %>%
        pivot_longer(-cit) %>%
        inner_join(uniq_dat, by=c('name' = 'idx')) %>%
        left_join(
             calc_post_strata(metadata, strata_cols) %>% filter(type == 'participant-visits'),
             by=strata_cols) %>%
         group_by(cit, comm_type) %>%
         summarise(
                value = sum(value*n)/sum(n), .groups='drop') %>%
         rename(name = comm_type)
    bounds = setNames(
        prop_mi %>% 
            group_by(name) %>% 
            group_map(~calc_hpd_hist_bounds(.x$value)),
        unlist(prop_mi %>%
            group_by(name) %>% 
            group_map(~.y$name)))
    # todo make multi shadded hist function
    # histogram bin width
    bw = clean_numeric((max(prop_mi$value) - min(prop_mi$value))/approx_bins)
    ll = 0
    ul = ceiling(max(prop_mi$value)/bw)*bw
    breaks = seq(ll, ul, bw)
    prop_mi$bin = cut(as_vector(prop_mi$value), breaks, dig.lab=10)
    prop_mi_sum = prop_mi %>%
        group_by(name, bin) %>%
        summarise(n=n(), .groups='drop') %>%
        group_by(name) %>%
        mutate(d = n / sum(n)) %>%
        ungroup() %>% 
        mutate(
                min = as.numeric(str_replace(str_split(bin, ",", simplify=TRUE)[,1], "\\(", "")),
                max = as.numeric(str_replace(str_split(bin, ",", simplify=TRUE)[,2], "\\]", "")),
                x =  (min + max)/2) %>%
        rowwise() %>%
        mutate(
                min_shade_bin = cut(min, bounds[[name]], labels=FALSE),
                max_shade_bin = cut(max, bounds[[name]], labels=FALSE)) %>%
        ungroup() %>%
        mutate(shade_bin = case_when(
                    # if they match, just choose one
                    min_shade_bin == max_shade_bin ~ min_shade_bin,
                    # if they don't match choose the one that pushes the interval out the most
                    (min_shade_bin != max_shade_bin) & 
                        min_shade_bin < 3 ~ max_shade_bin,
                    (min_shade_bin != max_shade_bin) & 
                        min_shade_bin >= 3 ~ min_shade_bin),
            shade_cat = paste(name, shade_bin, sep='_'))
    p = ggplot(prop_mi_sum, aes(x=x*100, y=d, fill=shade_cat, group=name))+
        geom_bar(color='#333333', stat="identity", position="identity") + 
        scale_fill_manual(values=colors_dict,
            breaks=c("inland_3", "fishing_3"),
            labels=c('inland', 'fishing'),
            name=NULL) +
        guides(
            fill = guide_legend(position="inside")) + 
        xlab(expression(atop('prevalence of', paste('multiple infections (', bar(delta), ', %)')))) +
        ylab('density') +
        xlim(0,NA)+
        scale_y_continuous(breaks=NULL, name='density', expand = expansion(mult = c(0, .15))) +
        gtheme +
        theme(legend.justification.inside=c(1,1), 
            axis.title=element_text(size=14))
    return(p)
}


plot_comm_type_sexpever = function(draws, dm, dat, colors_dict, plot_age_cat = '(24,29]'){
    # read in standardization values for each community
    age_cat_std = dat %>% filter(finalhiv == 'P' & sex == 'M' & age_cat_fine == plot_age_cat) %>%
        select(comm_type, plhiv_sexpever_mean) %>%
        unique()
    # plotting range
    sexpever = seq(1,30,0.5)
    # standardardized sexpever across plotting range by community type
    sexpever_std = list(
        inland = sexpever - 
            (age_cat_std %>% filter(comm_type == 'inland'))$plhiv_sexpever_mean,
        fishing = sexpever - 
            (age_cat_std %>% filter(comm_type == 'fishing'))$plhiv_sexpever_mean)
    # generate simulated design matrices
    # first two columns are for plotting not for DM
    sim_dm = tibble(
            comm_type = c(
                rep('fishing', length(sexpever_std$fishing)),
                rep('inland', length(sexpever_std$inland))),
            sexpever = c(sexpever, sexpever),
            intercept = 1,
            comm_typefishing = 
                c(rep(1, length(sexpever_std$fishing)),
                    rep(0, length(sexpever_std$inland))),
            comm_typeinland = 
                c(rep(0, length(sexpever_std$fishing)),
                    rep(1, length(sexpever_std$inland))),
            plhiv_sexpever_std = c(sexpever_std$fishing, sexpever_std$inland)) %>%
        mutate(
            `comm_typefishing*plhiv_sexpever_std` = comm_typefishing*plhiv_sexpever_std,
            `comm_typeinland*plhiv_sexpever_std` = comm_typeinland*plhiv_sexpever_std)
    # get coefficients 
    draws = draws %>% mutate(cit = seq(1,n()))
    coeffs = draws[c('cit', 'logit_prob_mi_baseline', colnames(draws)[grepl('logit_prob_mi_coeffs', colnames(draws))])]
    colnames(coeffs)[3:ncol(coeffs)] = dm$col
    # finally, calculate prob_mi
    prob_mi = tibble(bind_rows(as_tibble(apply(as.matrix(coeffs)[,2:ncol(coeffs)], 1, 
            function(x) inv_logit(colSums(t(sim_dm[3:ncol(sim_dm)])*x)))) %>%
        mutate(dm_idx = seq(1,n())) %>%
        pivot_longer(-dm_idx, names_to='cit') %>%
        left_join(sim_dm %>% mutate(
                dm_idx = seq(1,n())),
            by='dm_idx') %>%
        group_by(comm_type, sexpever, plhiv_sexpever_std) %>%
        group_map(~cbind(.y, 
            tibble(median = median(.x$value)),
            get_hpd(.x$value) %>%
                rename(
                    'lower95' = lower,
                    'upper95' = upper),
            get_hpd(.x$value, credMass=0.5) %>%
                rename(
                    'lower50' = lower,
                    'upper50' = upper)))))
    # finally, we can plot
    age_cat_split = str_split(gsub('\\(', '', gsub('\\]', '', plot_age_cat)), ',', simplify=TRUE)
    title = paste(c('men, ', as.numeric(age_cat_split[1])+1, ' to ', age_cat_split[2], ' years old'), collapse='')
    p = ggplot(prob_mi, aes(x=sexpever, y=median*100, ymin=lower95*100, ymax=upper95*100, color=comm_type, fill=comm_type)) + 
        geom_ribbon(alpha=0.25, linewidth=0.25) +
        geom_ribbon(aes(ymin=lower50*100, ymax=upper50*100), alpha=0.75, linewidth=0.25) +
        geom_line() +
        xlab('lifetime sex partners') +
        ylim(0, ceiling(max(prob_mi$upper95)*100+1)) + 
        ylab(expression(atop('predicted prevalence of', paste('multiple infections (', bar(delta), ', %)')))) +
        scale_fill_manual(values=colors_dict, guide=NULL) +
        scale_color_manual(values=
            c(
                inland = 
                    colors_dict[['inland_dark']],
                fishing = 
                    colors_dict[['fishing_dark']]),
            guide=NULL) +
        gtheme +
        theme(axis.title=element_text(size=14)) +
        ggtitle(title)
return(p)
}

p = arg_parser("plot summary of simulated data")
p = add_argument(p, "--dat", help="input data file")
p = add_argument(p, "--metadata", help="input metadata file")
p = add_argument(p, "--filter", help="data filter", default="id_subgraph_reads > 0 & window_type == 'unique'", nargs=1)
p = add_argument(p, "--commTypeFit", help="stan fit object", nargs=1)
p = add_argument(p, "--sexpeverFit", help="stan fit object", nargs=1)
p = add_argument(p, "--varSelectFit", help="stan fit object", nargs=1)
p = add_argument(p, "--out", help="output figure name", nargs=1)
p = add_argument(p, "--colors", help="tsv file with color codes", nargs=1, default='config/colors.tsv')
args = parse_args(p)

#args$commTypeFit = "fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_age_sex_comm.Rds"
#args$sexpeverFit = "fit/211220_allreads_phsc_all_subgraphs_format_par_m_deep-phyloMI_sexpever_men.Rds"
#args$varSelectFit = "fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_var_select.Rds"
#args$dat = "output/211220_allreads_phsc_all_subgraphs_format_par.tsv"
#args$out = "empirical_mi_predictors"
#args$metadata = 'output/211220_allreads_phsc_metadata.tsv'
# pre-defined colors
cols = read_tsv(args$colors, show_col_types=FALSE)
args$colors_dict = setNames(cols$color, cols$var)
colors_dict = args$colors_dict

# read in data
dat = summarise_n_d(tabulate_n_d(read_tsv(args$dat, show_col_types = FALSE) %>% 
    filter(eval(parse(text=args$filter))))) 
metadata = read_tsv(args$metadata, show_col_types=FALSE)

# read in and process community type fit
comm_type_fit = readRDS(args$commTypeFit)
comm_type_fit_draws = as_tibble(comm_type_fit$draws(
        inc_warmup = FALSE,
        format = "draws_df"))

p1 = plot_comm_type(comm_type_fit_draws, dat, metadata, args$colors_dict)

# read in and process sexpever fit
sexpever_fit = readRDS(args$sexpeverFit)
sexpever_fit_draws = as_tibble(sexpever_fit$draws(
        inc_warmup = FALSE,
        format = "draws_df"))

# read in design matrix
sexpever_dm = read_tsv(gsub('.Rds', '_design_cols.tsv', args$sexpeverFit), show_col_types=FALSE) %>%
    filter(type == 'mi')


p2 = plot_comm_type_sexpever(sexpever_fit_draws, sexpever_dm, dat, colors_dict)


## read in and process all fit
all_fit = readRDS(args$varSelectFit)

all_fit_draws = as_tibble(all_fit$draws(
        inc_warmup = FALSE,
        format = "draws_df"))

all_mi_dm = read_tsv(gsub('.Rds', '_design_cols.tsv', args$varSelectFit), show_col_types=FALSE) %>%
    filter(type == 'mi')

p3 = plot_coeffs(all_fit_draws,all_mi_dm, colors_dict)


p = plot_grid(
    plot_grid(p1, NULL, p2, nrow=1, rel_widths=c(0.5,0.025,0.5), 
        labels = c('a', '', 'b'), label_colour = "#333333"),
    NULL, p3, ncol=1, rel_heights=c(0.4, 0.025, 0.6), 
        labels=c('', '', 'c'), label_colour = "#333333")

ggsave(paste(c('figures/', args$out, '.pdf'), collapse=''), p, width=8, height=10)


