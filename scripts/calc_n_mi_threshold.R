suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(require(HDInterval))
suppressMessages(source('scripts/utils.R'))


calc_hpd_bounds = function(d){
		suppressMessages(require(HDInterval))
		bounds95 = hdi(d, credMass=0.95)
		bounds50 = hdi(d, credMass=0.50)
		return(c(low95 = bounds95[[1]], low50 = bounds50[[1]], 
			median=median(d), high50=bounds50[[2]], high95=bounds95[[2]]))
	}

# for each cutoff, get median 95% HPD and 50% HPD 
get_n_mi_hpd = function(p_dat, cutoff){
	return(c(calc_hpd_bounds((p_dat %>%
			mutate(mi = value >= cutoff) %>%
			group_by(cit) %>%
			summarise(n_mi = sum(mi)))$n_mi), cutoff=cutoff))
}

#### MAIN CODE BLOCK ####
p <- arg_parser("calc mi dichotomization thresholds")
p <- add_argument(p, "--fit", help="stan fit object", nargs=1)
p <- add_argument(p, "--dat", help="input data file", nargs=1)
p <- add_argument(p, "--out", help="output figure name", nargs=1)
args <- parse_args(p)

#args$fit = 'fit/211220_allreads_phsc_all_subgraphs_format_par_full_model_.Rds'
#args$dat = 'output/211220_allreads_phsc_all_subgraphs_format_par.tsv'

dat = summarise_n_d(tabulate_n_d(read_tsv(args$dat, show_col_types = FALSE) %>% 
	filter(window_type == 'unique'))) 
fit = readRDS(args$fit)

#  number of detected MIs by posterior cut-off
cutoffs = seq(0.01,1,0.01)
param=c('prob_MI', 'ind_log_prob_mi')

fit_draws = as_tibble(fit$draws(
		variables = c(param),
        inc_warmup = FALSE,
        format = "draws_df"))

fit_draws = fit_draws[,(!grepl('ppr', colnames(fit_draws)) & 
    !grepl('prob_seq_any', colnames(fit_draws)) & 
    !grepl('prob_seq_1', colnames(fit_draws)) & 
    !grepl('prob_seq_MI', colnames(fit_draws)) & 
    !grepl('N_nnzero_times_log_prob_MI', colnames(fit_draws)) & 
    !grepl('X_mi_missing', colnames(fit_draws)))] %>%
  mutate(cit = seq(1,n())) %>%
  pivot_longer(-cit)

p_dat = bind_rows(
	fit_draws %>%
		mutate(value = exp(value)) %>% 
	    group_by(name) %>% 
	    group_map(~get_hpd(.x$value) %>% 
	    mutate(name=.y$name))) %>%
	rename(`0.025` = lower, `0.975` =upper) %>%
	left_join(
		bind_rows(
			fit_draws %>%
				mutate(value = exp(value)) %>% 
			    group_by(name) %>% 
			    group_map(~get_hpd(.x$value, credMass=0.50) %>% 
			    mutate(name=.y$name))) %>%
		rename(`0.25` = lower, `0.75`=upper),
		by='name') %>%
	left_join(
		fit_draws %>%
			mutate(value = exp(value)) %>% 
		    group_by(name) %>%
		    summarise(`0.50`=median(value)),
		by='name') %>%
	filter(name != '.chain' & name != '.draw' & name != '.iteration') %>%
	mutate(name = gsub('_log', '', name))

p_dat = p_dat %>% 
	rename(id = name) %>%
	pivot_longer(-id)

cutoffs = seq(0.01, 1, 0.001)
n_mi_cutoff = bind_rows(lapply(cutoffs, 
	function(x){p_dat %>% group_by(name) %>% 
		summarise(n_mi = sum(value >= x)) %>% mutate(cutoff = x)}))

n_mi_cutoff = n_mi_cutoff %>% pivot_wider(names_from=name, values_from=n_mi)

# estimated median prevalence of multiple infections
estimated_n_mi = (fit_draws %>% filter(name == 'prob_MI') %>% summarise(m = median(value)))$m * 
	length(unique(p_dat$id))

n_mi_threshold = n_mi_cutoff %>% pivot_longer(-cutoff) %>%
	filter(value <= estimated_n_mi) %>%
	group_by(name) %>% filter(cutoff == min(cutoff)) 


out = bind_rows(
	n_mi_cutoff %>%
		mutate(cutoff = as.character(cutoff)), 
	n_mi_threshold %>%
		mutate(cutoff = 'n_mi') %>%
		pivot_wider(names_from=name, values_from=value) %>%
		mutate(cutoff = as.character(cutoff)),
	n_mi_threshold %>%
		select(-value) %>%
		rename(value=cutoff) %>%
		mutate(cutoff = 'empirical') %>%
		pivot_wider(names_from=name, values_from=value) %>%
		mutate(cutoff = as.character(cutoff)))


write_tsv(out, str_replace(args$fit, '.Rds', '_n_mi_dichotomized.tsv'))

