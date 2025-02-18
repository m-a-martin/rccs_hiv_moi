suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(require(cowplot))
suppressMessages(require(patchwork))
suppressMessages(require(HDInterval))
suppressMessages(source('scripts/utils.R'))


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
p = add_argument(p, "--sexpeverFit", help="stan fit object", nargs=1)
p = add_argument(p, "--out", help="output figure name", nargs=1)
p = add_argument(p, "--colors", help="tsv file with color codes", nargs=1, default='config/colors.tsv')
args = parse_args(p)

#args$sexpeverFit = 'fit/211220_allreads_phsc_all_subgraphs_format_par_m_extended_model_sexpever_men.Rds'
#args$dat = 'output/211220_allreads_phsc_all_subgraphs_format_par.tsv'

cols = read_tsv(args$colors, show_col_types=FALSE)
args$colors_dict = setNames(cols$color, cols$var)
colors_dict = args$colors_dict

# read in and process sexpever fit
sexpever_fit = readRDS(args$sexpeverFit)
sexpever_fit_draws = as_tibble(sexpever_fit$draws(
        inc_warmup = FALSE,
        format = "draws_df"))

# read in design matrix
sexpever_dm = read_tsv(gsub('.Rds', '_design_cols.tsv', args$sexpeverFit), show_col_types=FALSE) %>%
	filter(type == 'mi')

# we need to calculate prop_mi in inland and fishing communities
# as a function of standardized sexpever 
# we need standardization values
dat = read_tsv(args$dat, show_col_types = FALSE)


plot_list = list()
for (age_cat in sort(unique(dat$age_cat_fine))){
	plot_list[[length(plot_list)+1]] = 
		plot_comm_type_sexpever(sexpever_fit_draws, sexpever_dm, dat, 
			colors_dict, plot_age_cat = age_cat) +
		ylim(0, 32) +
		ylab(NULL) 
}

plot_list[[1]] = plot_list[[1]] + 
	ylab(expression(atop('predicted prevalence of', paste('multiple infections (', bar(delta), ', %)')))) +
	guides(color = guide_legend(position="inside", title=NULL),
		fill = guide_legend(position="inside", title=NULL)) +
	theme(legend.justification.inside = c(0, 1))

p = plot_grid(plotlist=plot_list, nrow=1, rel_widths=c(1,c(rep(0.8, length(plot_list)-1))))

ggsave(paste(c('figures/', args$out, '.pdf'), collapse=''), p, width=22, height=4.8)

