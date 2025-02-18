suppressMessages(require(cmdstanr))
suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(require(HDInterval))
suppressMessages(require(posterior))
suppressMessages(source('scripts/utils.R'))

summarize_group = function(x, name){
  return(get_hpd(x$value) %>%
              mutate(
                name = name,
                median = median(x$value),
                bulk_ess = ess_bulk(x$value),
                tail_ess = ess_tail(x$value),
                rhat = rhat(x$value)))
}


get_coeffs = function(fit_draws){
  baseline = fit_draws[c('cit', 'logit_prob_mi_baseline')]
  # have to account for shrinkage priors when exist
  if (any(grepl('w\\[1\\]', colnames(fit_draws)))){
    coeffs = fit_draws[grepl('w\\[', colnames(fit_draws))]
  }else{
    coeffs = fit_draws[grepl('logit_prob_mi_coeffs\\[', colnames(fit_draws))]
  }
  coeffs = bind_cols(baseline, coeffs)
  return(coeffs)
}


#### --------------------- ####
#### MAIN CODE BEGINS HERE ####
#### --------------------- ####
p <- arg_parser("calculate post stratification. assumes no continuous variables in design matrix!")
#### METADATA ####
p <- add_argument(p, "--metadata", help="input data file", nargs=1)
#### STAN FIT ####
p <- add_argument(p, "--fit", help="stan fit object", nargs=1)
#### STRATIFICATION VARS ####
p = add_argument(p, "--strataVars", nargs=Inf, 
  help='variables for design matrix and post stratification. assumes same for both, assumes matches design matrix for model fit.')
#### RISKS AND RISKS RATIOS ####
p = add_argument(p, "--strataPrevs", nargs=Inf, 
  help='post hoc calculation of risk based on design matrix. cannot be more detailed than poststrata risks')
p = add_argument(p, "--strataPrevRatios", nargs=Inf,  
  help='post hoc calculation of risk ratios based on design matrix')
args <- parse_args(p)
#args$metadata = 'output/211220_allreads_phsc_metadata.tsv'
#args$fit = 'fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_age_sex_comm.Rds'
#args$strataVars = c("age_cat_coarse", "sex", "comm_type")
#args$strataPrevs = c("TRUE", "comm_type == 'fishing'", "comm_type == 'inland'")
#args$strataPrevRatios = c("comm_type == 'fishing':comm_type=='inland'")

# 1. CALC POSTSTRATIFICATION #
metadata = read_tsv(args$metadata, show_col_types=FALSE)
# get post-strata N
post_strata = calc_post_strata(metadata, args$strataVars) %>%
	select(-p) %>%
	filter(type == "participant-visits") %>%
	select(-type)

# 2. GENERATE DESIGN MATRIX # 
uniq_metadata = metadata %>% 
	filter(phsc_par == TRUE) %>%
	select(args$strataVars) %>% 
	unique() %>%
	left_join(post_strata, by=args$strataVars) %>%
	mutate(idx = seq(1,n()))

if ('round' %in% colnames(uniq_metadata)){
	metadata = metadata %>% mutate(round = as.character(round))
}

dm = generate_dm(uniq_metadata, args$strataVars)$dm

# 3. PROCESS FIT AND PREPARE TO CALCULATE OUTCOME #
fit = readRDS(args$fit)

fit_draws = as_tibble(fit$draws(
				variables = c('logit_prob_mi_baseline', 'logit_prob_mi_coeffs'),
        inc_warmup = FALSE,
        format = "draws_df")) %>%
  mutate(cit = seq(1,n()))

# get coefficients
coeffs = get_coeffs(fit_draws)
# for each row in unique data frame calculate prob_mi
prob_mi = apply(as.matrix(coeffs)[,2:ncol(coeffs)], 1, 
        function(x) inv_logit(colSums(t(cbind(1, dm))*x)))
colnames(prob_mi) = as.character(seq(1,ncol(prob_mi)))
prob_mi =  as_tibble(prob_mi) %>% 
	mutate(idx = seq(1,n())) %>%
	pivot_longer(-idx) %>%
	mutate(cit = as.numeric(name)) %>%
	select(-name) %>%
	left_join(uniq_metadata, by='idx')

saveRDS(prob_mi, file = str_replace(args$fit, ".Rds", "_poststrat.Rds"))

# 4. CALCULATE DESIRED OUTCOME MEASURES #
poststrat_sum = prob_mi %>% 
	group_by(cit) %>%
	summarise(value = sum(value*n)/sum(n)) %>%
	summarise(summarize_group(., 'wprob_mi'))

for (i in args$strataPrevs){
	# get indices corresponding to strata
	strata_prob_mi = uniq_metadata %>%
		filter(eval(parse(text=i)))
	poststrat_sum = bind_rows(
		poststrat_sum,
		prob_mi %>% filter(eval(parse(text=i))) %>%
			group_by(cit) %>%
			summarise(value = sum(value*n)/sum(n)) %>%
			summarise(summarize_group(., paste(c('wprob_mi_', gsub('"', "", i)), collapse=''))))
}


for (i in args$strataPrevRatios){
	i_split = str_split(i, ":", simplify=TRUE)
	poststrat_sum = bind_rows(
		poststrat_sum,
		prob_mi %>% 
	        filter(eval(parse(text=i_split[1]))) %>%
	        group_by(cit) %>%
	        summarise(num = sum(value*n)/sum(n)) %>%
	        left_join(
	          prob_mi %>% 
	            filter(eval(parse(text=i_split[2]))) %>%
	            group_by(cit) %>%
	            summarise(denom = sum(value*n)/sum(n)),
	          by='cit') %>%
	        mutate(value = num/denom) %>%
	        summarise(summarize_group(., paste(c('wpr_mi', "_", gsub('"', "", i)), collapse=''))))
}


write_tsv(poststrat_sum, str_replace(args$fit, ".Rds", "_poststrat_summary.tsv"))

