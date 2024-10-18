suppressMessages(library(tidyverse))
suppressMessages(library(argparser))
suppressMessages(source('scripts/utils.R'))


simulate_data = function(args){
	set.seed(42L)
	# set up basic data
	sd = tibble(	
			id = 1:args$N_ind,
			log_VL_obs = rnorm(args$N_ind, 
				args$log_VL_obs_mean, args$log_VL_obs_sd)) %>%
			mutate(
				scaled_log10_vl_obs = c(scale(log_VL_obs)),
				seq_success = inv_logit(args$logit_prob_seq_baseline + 
					args$logit_seq_prob_lvl_effect*scaled_log10_vl_obs + 
					rnorm(args$N_ind, 0, args$logit_prob_seq_ind_sd)))
	# initialize multiple infection status
	tmp_mi = rep(0, nrow(sd))
	# add risk factor status
	if (args$prob_risk_factor > 0){
		# assign risk factor
		# ensures prob_risk factor is exact
		tmp = rep(0, nrow(sd))
		tmp[sample(1:nrow(sd), 
      		floor(nrow(sd) * args$prob_risk_factor), 
      		replace=F)] = 1
		sd[paste("risk_factor_MI_", 1, sep='')] = tmp
		# how many individuals have risk factor
		n_w_risk = sum(tmp)
		# how many individuals do not have risk factor
		n_wo_risk = sum(tmp == 0)
		tmp_mi[tmp == 1][sample(1:n_w_risk,
			floor(args$prob_MI_risk_factor*n_w_risk), 
			replace=F)] = 1
		tmp_mi[tmp == 0][sample(1:n_wo_risk,
			floor(args$prob_MI_baseline*n_wo_risk), 
			replace=F)] = 1
		# if null risk factors then add those
		if (args$null_risk_factors > 0){
			for (i in 1:args$null_risk_factors){
				tmp_null = rep(0, nrow(sd))
				tmp_null[sample(1:nrow(sd), 
		      		floor(nrow(sd) * args$prob_risk_factor), 
		      		replace=F)] = 1
				sd[paste("risk_factor_MI_", (i + as.integer(args$prob_risk_factor > 0)), sep='')] = tmp_null
			}	
		}
	}else{
		tmp_mi[sample(1:nrow(sd), floor(nrow(sd) * args$prob_MI_baseline), replace=F)] = 1
	}
	sd["has_MI"] = tmp_mi
	simulate_window_dat = function(N, seq_success, fpr, fnr, has_MI){
		v1 = as.integer(runif(N) < seq_success)
		v2 = rep(0, N)
		fp1 = as.integer(runif(N) < fpr)
		fp2 = rep(0, N)
		if (has_MI == 1){
			v2 = as.integer(runif(N) < seq_success)
			fp2 = as.integer(runif(N) < fpr)
		}
	 	fn <- as.integer(runif(N) < fnr) * v1 * v2
	 	return(c(
	 			"N_obs"=sum(v1 | v2),
	 			"MI_obs"=sum((v1 * v2 * (1-fn)) | (v1*(1-v2)*fp1) | ((1-v1)*v2*fp2))))
	 	}
	sd$N_windows = args$N_windows
	sd$N_obs = 0
	sd$MI_obs = 0
	for (i in 1:args$N_ind){
		N_MI = simulate_window_dat(args$N_windows, sd$seq_success[i], args$prob_MI_fpr, args$prob_MI_fnr, sd$has_MI[i])	
		sd$N_obs[i] = N_MI["N_obs"]
		sd$MI_obs[i] = N_MI["MI_obs"]
	}
	# remove samples with no sequenced data
	sd = sd %>% filter(N_obs > 0)
	# renumber individuals
	sd = sd %>% mutate(id = seq(1,n()))
	return(sd)
}


p <- arg_parser("simulate data")
p <- add_argument(p, "--N_ind", help="number of individuals", 
	nargs=1, type='integer', default=2000)
p <- add_argument(p, "--N_windows", help="number of genome windows", 
	nargs=1, type='integer', default=29)
p <- add_argument(p, "--logit_prob_seq_baseline", help="baseline log(odds) of sequencing success", 
	nargs=1, type='double', default=2)
p <- add_argument(p, "--logit_seq_prob_lvl_effect", 
	help="increase in log(odds) of sequencing success per unit increase in normalized viral load", 
	nargs=1, type='double', default=2)
p <- add_argument(p, "--logit_prob_seq_ind_sd", 
	help="standard deviation on inidividual variation of sequencing success", 
	nargs=1, type='double', default=1)
p <- add_argument(p, "--prob_MI_baseline", help="population prevalence of multiple infection", 
	nargs=1, type='double', default=0.05)
p <- add_argument(p, "--prob_MI_fpr", help="multiple infection false positive rate", 
	nargs=1, type='double', default=0)
p <- add_argument(p, "--prob_MI_fnr", help="multiple infection false negative rate", 
	nargs=1, type='double', default=0)
p <- add_argument(p, "--prob_risk_factor", help="population prevalence of harboring MI risk factor", 
	nargs=1, type='double', default=0)
p <- add_argument(p, "--prob_MI_risk_factor", help="population prevalence of multiple infections among those with risk factor", 
	nargs=1, type='double', default=0)
p <- add_argument(p, "--null_risk_factors", help="number of risk factors unassociated with multiple infection status", 
	nargs=1, type='integer', default=0)
p <- add_argument(p, "--log_VL_obs_mean", help="mean log10 viral load", 
	nargs=1, type='float', default=1e3)
p <- add_argument(p, "--log_VL_obs_sd", help="sd log10 viral load", 
	nargs=1, type='float', default=1e5)
p <- add_argument(p, "--out", help="output file name", 
	nargs=1, type='character')
args <- parse_args(p)

args$logit_prob_MI = logit(args$prob_MI_baseline)

if (!file.exists('simulations/')){
	dir.create('simulations')	
}

sd = simulate_data(args)
op = paste(c('simulations/', args$out, '.tsv'), collapse='')
cat(op)

op_base = str_split(op, '/', simplify=TRUE)
op_base = paste(op_base[1:length(op_base)-1], collapse='/')
if (!dir.exists('simulations')){dir.create('simulations')}
if (!dir.exists(op_base)){dir.create(op_base)}

write_tsv(sd, op)

out_arg = names(args)[4:length(args)]
out_val = unlist(args)[4:length(args)]

out_arg = c(out_arg, "logit_prob_seq_coeffs[1]")
out_val = c(out_val, args$logit_seq_prob_lvl_effect)

if (args$prob_risk_factor > 0){
	out_arg = c(out_arg ,'w[1]')
	# calculate odds ratio
	out_val = c(out_val, logit(args$prob_MI_risk_factor) - logit(args$prob_MI_baseline))
	for (i in 1:args$null_risk_factors){
		out_arg = c(out_arg, paste(c('w[', as.character(i+1), ']'), collapse=''))
		out_val = c(out_val, 0)
	}
}

if (!dir.exists('simulations')){
	dir.create('simulations')
}

write_tsv(tibble(arg = out_arg, value=out_val),
	paste(c('simulations/', args$out, '_params.tsv'), collapse=''))
