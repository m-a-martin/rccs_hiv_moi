suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(require(haven))
suppressMessages(require(HDInterval))
suppressMessages(source('scripts/utils.R'))

nums = c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten")
calc_basic_stats = function(out, metadata, round_digits=2){
	# unpack input
	labels = out[["labels"]]
	vals = out[['vals']]
	#### BASIC STATISTICS ABOUT DATA SET ####
	# number of survey rounds
	labels = c(labels, "n_rounds")
	vals = c(vals, length(unique(metadata$round)))
	labels = c(labels, "n_rounds_char")
	vals = c(vals, nums[length(unique(metadata$round))])
	labels = c(labels, "min_round")
	vals = c(vals, min(metadata$round))
	labels = c(labels, "max_round")
	vals = c(vals, max(metadata$round))
	# number of survey participants 
	labels = c(labels, "n_participants")
	vals = c(vals, format(nrow(metadata %>% select(study_id) %>% unique()), big.mark=','))
	n_participant_visits = nrow(metadata %>% select(study_id, round) %>% unique())
	labels = c(labels, "n_participant_visits")
	vals = c(vals, format(n_participant_visits,big.mark=','))
	# HIV+ participants
	labels = c(labels, "n_hiv_participants")
	vals = c(vals, format(nrow(metadata %>% filter(finalhiv == "P") %>% 
		select(study_id) %>% unique()), big.mark=','))
	labels = c(labels, "p_hiv_participants")
	vals = c(vals, round(100*nrow(metadata %>% filter(finalhiv == "P") %>% 
		select(study_id) %>% unique()) / nrow(metadata %>% 
		select(study_id) %>% unique()),round_digits))
	# HIV+ participant visits
	n_hiv_participant_visits = 
		nrow(metadata %>% filter(finalhiv == "P"))
	labels = c(labels, "n_hiv_participant_visits")
	vals = c(vals, format(n_hiv_participant_visits, big.mark=','))
	labels = c(labels, "p_hiv_participant_visits")
	vals = c(vals, round(100*n_hiv_participant_visits / n_participant_visits, round_digits))
	labels = c(labels, "n_neg_participant_visits")
	vals = c(vals, format(nrow(metadata %>% filter(finalhiv != "P")), big.mark=','))
	# Viremic participants
	labels = c(labels, "n_viremic_participants")
	vals = c(vals, format(nrow(metadata %>% filter(finalhiv == "P" &  log10_copies >= 3) %>% 
		select(study_id) %>% unique()), big.mark=','))
	# PHSC participants
	labels = c(labels, "n_phsc_participants")
	vals = c(vals, format(nrow(metadata %>% filter(finalhiv == "P" &  log10_copies >= 3 & phsc==TRUE) %>% 
		select(study_id) %>% unique()), big.mark=','))
	labels=c(labels, 'min_survey_date')
	vals = c(vals, as.character(metadata$min_survey_date[[1]]))
	labels=c(labels, 'max_survey_date')
	vals = c(vals, as.character(metadata$max_survey_date[[1]]))
	labels=c(labels, 'min_survey_year')
	vals = c(vals, str_split(metadata$min_survey_date[[1]], " ", simplify=T)[,2])
	labels=c(labels, 'max_survey_year')
	vals = c(vals, str_split(metadata$max_survey_date[[1]], " ", simplify=T)[,2])
	# dates per round
	survey_round_years = metadata %>%
		select(round, round_median_year) %>%
		unique()
	for (i in 1:nrow(survey_round_years)){
		labels = c(labels, paste(c('round_', survey_round_years$round[i], '_median_year'), collapse=''))
		vals = c(vals, survey_round_years$round_median_year[i])
	}
	# repack output
	out[["labels"]] = labels
	out[["vals"]] = vals
	return(out)
}


calc_epi_stats = function(out, metadata){
	format_p_value = function(x){
		if (x < 1E-4){
			return("<1.0\\\\times10^{-4}")
		}else{
			return(sprintf("%.4f", round(x,4)))
		}
	}
	# subset to just participants from participant-visits
	# if included in phsc, then use that visit
	phsc_pars = metadata %>% filter(phsc_par == TRUE)
	non_phsc_pars = metadata %>% 
		group_by(study_id) %>%
		# here get only people not included in phyloscanner analysis
		filter(all(phsc_par == FALSE | is.na(phsc_par))) %>%
		# get highest viral load
		filter(log10_copies == max(log10_copies) | is.na(log10_copies)) %>%
		# if ties get first
		filter(round == min(round)) %>%
		ungroup()
	par_metadata = bind_rows(phsc_pars, non_phsc_pars)
	# unpack input
	labels = out[["labels"]]
	vals = out[['vals']]
	#### EPI DATA (TABLE 1) ####
	epi_dat = list()
	epi_dat[['participant']] = par_metadata
	epi_dat[['hiv_participant']] = par_metadata %>% 
		filter(finalhiv == "P")
	epi_dat[['viremic_participant']] = par_metadata %>% 
		filter(finalhiv == "P" & log10_copies >= 3)
	epi_dat[['seq_viremic_participant']] = metadata %>% 
		filter(finalhiv == "P" & log10_copies >= 3 & phsc_par == TRUE)
	for (var in c('round', 'sex', 'age_cat_coarse', 'comm_type', 'log10_copies_cat')){
		tabulated = tibble()
		for (dat in names(epi_dat)){
			sum_dat = epi_dat[[dat]] %>% 
				group_by_at(var) %>% summarise(n=n(), .groups="drop") %>%
				mutate(p = 100*(n / sum(n)) )
			tabulated = bind_rows(tabulated, sum_dat %>% mutate(dat=dat))
			n_lab=paste(c('n', dat, var), collapse='_')
			p_lab=paste(c('p', dat, var), collapse='_')
			for (idx in seq(1,nrow(sum_dat))){
				labels = c(labels, paste(n_lab, as_vector(sum_dat[,1])[idx], sep='_'))
				vals = c(vals, format(as_vector(sum_dat[,2])[idx], big.mark=','))
				labels = c(labels, paste(p_lab, as_vector(sum_dat[,1])[idx], sep='_'))
				vals = c(vals, round(as_vector(sum_dat[,3])[idx], round_digits))
				
			}
		}
		# all v. hiv proportion and CI
		if (var != 'log10_copies_cat'){
			all_hiv_dat = tabulated %>% select(-p) %>% 
				filter(dat == 'participant' | dat == 'hiv_participant') %>%
				pivot_wider(names_from=dat, values_from=n) %>%
				mutate(
					p = hiv_participant / participant,
					n_tilde = participant + qnorm(0.025)^2,
					p_tilde = (1/n_tilde)*(hiv_participant + qnorm(0.025)^2/2),
					low = p_tilde + qnorm(0.025)*sqrt((p_tilde/n_tilde) * (1-p_tilde)),
					high = p_tilde - qnorm(0.025)*sqrt((p_tilde/n_tilde) * (1-p_tilde))) %>%
				select(-n_tilde, -p_tilde) %>%
				mutate(across(p:high, ~paste(round(.x*100, round_digits), "\\\\%", sep=''))) %>%
				unite("p_format", p:low, sep=' (') %>%
				unite("p_format", p_format:high, sep=' - ') %>%
				mutate(p_format = paste(p_format, ')', sep=''))
			labels = c(labels, paste(paste(c('p_sample_hiv_participant_', var, '_'), collapse=''),
				as.vector(all_hiv_dat[,1])[[1]], sep=''))
			vals = c(vals, as.vector(all_hiv_dat[,4])[[1]])
		}
		# viremic hiv v. sequenced viremic proportion and CI
		viremic_seq_dat = tabulated %>% select(-p) %>% 
			filter(dat == 'viremic_participant' | dat == 'seq_viremic_participant') %>%
			pivot_wider(names_from=dat, values_from=n) %>%
			mutate(
				p = seq_viremic_participant / viremic_participant,
				n_tilde = viremic_participant + qnorm(0.025)^2,
				p_tilde = (1/n_tilde)*(seq_viremic_participant + qnorm(0.025)^2/2),
				low = p_tilde + qnorm(0.025)*sqrt((p_tilde/n_tilde) * (1-p_tilde)),
				high = p_tilde - qnorm(0.025)*sqrt((p_tilde/n_tilde) * (1-p_tilde))) %>%
			select(-n_tilde, -p_tilde) %>%
			mutate(across(p:high, ~paste(round(.x*100, round_digits), "\\\\%", sep=''))) %>%
			unite("p_format", p:low, sep=' (') %>%
			unite("p_format", p_format:high, sep=' - ') %>%
			mutate(p_format = paste(p_format, ')', sep=''))
		labels = c(labels, paste(paste(c('p_sample_seq_participant_', var, '_'), collapse=''),
			as.vector(viremic_seq_dat[,1])[[1]], sep=''))
		vals = c(vals, as.vector(viremic_seq_dat[,4])[[1]])
	}
	# repack output
	out[["labels"]] = labels
	out[["vals"]] = vals
	return(out)
}


calc_supp_epi_stats = function(out, metadata, round_digits=2){
	# subset to just participants from participant-visits
	# if included in phsc, then use that visit
	phsc_pars = metadata %>% filter(phsc_par == TRUE)
	non_phsc_pars = metadata %>% 
		group_by(study_id) %>%
		# here get only people not included in phyloscanner analysis
		filter(all(phsc_par == FALSE | is.na(phsc_par))) %>%
		# get highest viral load
		filter(log10_copies == max(log10_copies) | is.na(log10_copies)) %>%
		# if ties get first
		filter(round == min(round)) %>%
		ungroup()
	par_metadata = bind_rows(phsc_pars, non_phsc_pars)
	# unpack input
	labels = out[["labels"]]
	vals = out[['vals']]
	##### MISSING VALUES #####
	n_participants = nrow(par_metadata)
	n_seq_viremic = nrow(par_metadata %>% 
		filter(finalhiv == "P" & log10_copies >= 3 & phsc_par == TRUE))
	vars = c('finalhiv_orig', 'sex', 'age_cat_fine', 'male_circumcision', 
			'comm_type', 'married', 'sexpever', 'in_migrant', 'barworker')
	stopifnot(all(vars %in% colnames(par_metadata)))
	for (var in vars){
		missing = sum(is.na(par_metadata[[var]]))
		labels = c(labels, paste(c('n_participant_', var, '_missing'), collapse=''))
		vals = c(vals, missing)
		labels = c(labels, paste(c('p_participant_', var, '_missing'), collapse=''))
		vals = c(vals, round(100*missing / n_participants, round_digits))
		seq_viremic_missing = sum(is.na((par_metadata  %>% 
			filter(finalhiv == "P" & log10_copies >= 3 & phsc == TRUE))[[var]]))
		labels = c(labels, paste(c('n_seq_viremic_participant_', var, '_missing'), collapse=''))
		vals = c(vals, seq_viremic_missing)
		labels = c(labels, paste(c('p_seq_viremic_participant_', var, '_missing'), collapse=''))
		vals = c(vals, round(100*seq_viremic_missing / n_seq_viremic, round_digits))
	}
	labels = c(labels, 'n_participant_finalhiv_orig_available')
	vals = c(vals, format(sum(!is.na(par_metadata[['finalhiv_orig']])), big.mark=','))
	labels = c(labels, 'n_participant_sexpever_categorical')
	vals = c(vals, format(nrow(par_metadata %>% filter(sexpever == 92 | sexpever == 93)), big.mark=','))
	labels = c(labels, 'p_participant_sexpever_categorical')
	vals = c(vals, round(100*nrow(par_metadata %>% filter(sexpever == 92 | sexpever == 93)) /
		nrow(par_metadata), round_digits))
	# median viral load
	#labels = c(labels, 'median_viremic_log10_copies')
	#median_viremic_log10_copies = (metadata %>% filter(finalhiv == 'P' & log10_copies >= 3) %>%
	#		summarise(m = log10(median(10**log10_copies))))[[1]]
	#vals = c(vals, 
	#	round(median_viremic_log10_copies,round_digits))
	#labels=c(labels, '10pct_median_viremic_log10_copies')
	#vals = c(vals, 
	#	round(log10(0.10*(10**median_viremic_log10_copies)), round_digits))
	# repack output
	out[["labels"]] = labels
	out[["vals"]] = vals
	return(out)
}



calc_vl_stats = function(out, round_rccs_dat, round_digits=2){
	# unpack input
	labels = out[["labels"]]
	vals = out[['vals']]
	n_hiv_participant_visits = 
		nrow(round_rccs_dat %>% filter(finalhiv == "P"))
	#### VIRAL LOAD DATA ####
	labels = c(labels, "n_copies_participant_visits")
	n_copies_participant_visits = 
		nrow(round_rccs_dat %>% filter(finalhiv == "P" & !is.na(log10_copies)))
	vals = c(vals, format(n_copies_participant_visits, big.mark=','))
	labels = c(labels, "p_copies_participant_visits")
	vals = c(vals, 
		round(100*n_copies_participant_visits / 
			n_hiv_participant_visits, round_digits))
	labels = c(labels, "n_viremic_participant_visits")
	n_viremic_participant_visits = 
		nrow(round_rccs_dat %>% filter(finalhiv == "P" & !is.na(log10_copies) & log10_copies >= 3))
	vals = c(vals, format(n_viremic_participant_visits, big.mark=','))
	labels = c(labels, "p_viremic_participant_visits")
	vals = c(vals, round(100*n_viremic_participant_visits/n_copies_participant_visits,round_digits))
	# missing VL in R14 through R15
	p_r14_r15_nonmissing_vl = round_rccs_dat %>% filter(finalhiv == 'P' & (round ==14 | round == 15)) %>%
		summarise(p = 100*sum(log10_copies_cat != 'missing')/n())
	labels = c(labels, 'p_r14_r15_nonmissing_vl')
	vals = c(vals, round(p_r14_r15_nonmissing_vl$p, round_digits))
	# missing VL in R16 through R19
	p_r16_r19_nonmissing_vl = round_rccs_dat %>% filter(finalhiv == 'P' & (round >= 16 & round <= 19)) %>%
		summarise(p = 100*sum(log10_copies_cat != 'missing')/n())
	labels = c(labels, 'p_r16_r19_nonmissing_vl')
	vals = c(vals, round(p_r16_r19_nonmissing_vl$p, round_digits))
	# repack output
	out[["labels"]] = labels
	out[["vals"]] = vals
	return(out)
}


calc_phsc_config = function(out, phsc_dat){
	# unpack input
	labels = out[["labels"]]
	vals = out[['vals']]
	#### PHYLOSCANNER CONFIG ####
	labels = c(labels, "phsc_window_size")
	vals = c(vals, phsc_dat$window_end[1] - phsc_dat$window_start[1] + 1)
	window_starts = sort(unique(phsc_dat$window_start))
	labels = c(labels, "phsc_step_size")
	vals = c(vals, window_starts[2] - window_starts[1])
	labels = c(labels, "n_all_windows")
	vals = c(vals, length(unique(phsc_dat$window_start)))
	labels = c(labels, "n_unique_windows")
	vals = c(vals, length(unique((phsc_dat %>% filter(window_type == "unique"))$window_start)))
	# repack output
	out[["labels"]] = labels
	out[["vals"]] = vals
	return(out)
}



calc_phsc_input = function(out, metadata){
	# unpack input
	labels = out[["labels"]]
	vals = out[['vals']]
	phsc_samples = metadata %>% filter(log10_copies >= 3 & phsc == TRUE)
	#### PHYLOSCANNER INPUT ####
	labels=c(labels, "n_phsc_samples")
	n_phsc_samples = nrow(phsc_samples %>% select(study_id, round) %>% unique())
	vals = c(vals, format(n_phsc_samples, big.mark=','))
	labels=c(labels, "n_nophsc_viremic_samples")
	vals = c(vals, format(nrow(metadata %>% 
		filter(phsc == FALSE & finalhiv == "P" & log10_copies >= 3) %>%
		select(study_id, round) %>% unique()), big.mark=','))
	labels=c(labels, "n_phsc_participants")
	vals = c(vals, format(nrow(phsc_samples %>% select(study_id) %>% unique()), big.mark=','))
	labels = c(labels, "min_phsc_date")
	vals = c(vals, phsc_samples$min_phsc_date[[1]])
	labels = c(labels, "max_phsc_date")
	vals = c(vals, phsc_samples$max_phsc_date[[1]])
	labels = c(labels, "min_phsc_year")
	vals = c(vals, str_split(phsc_samples$min_phsc_date[[1]], " ", simplify=TRUE)[2])
	labels = c(labels, "max_phsc_year")
	vals = c(vals, str_split(phsc_samples$max_phsc_date[[1]], " ", simplify=TRUE)[2])


	# group by sequencing technology and get counts
	seq_tech_counts = metadata %>% 
		filter(log10_copies > 3 & phsc_par == TRUE) %>%
		group_by(round, sequencing_technology) %>%
		summarise(n=n(), .groups='drop')
	for (i in seq(1,nrow(seq_tech_counts))){
		labels=c(labels, 
			paste(c('n_seq_', 
				seq_tech_counts[i,]$sequencing_technology, 
				'_round_', 
				seq_tech_counts[i,]$round), collapse=''))
		vals = c(vals, seq_tech_counts[i,]$n)
	}
	# repack output
	out[["labels"]] = labels
	out[["vals"]] = vals
	return(out)
}


calc_weight_stats = function(out, metadata, round_digits){
	labels = out[["labels"]]
	vals = out[['vals']]
	ps = calc_post_strata(metadata, c("age_cat_coarse", "sex", "comm_type"))
	ps_dat = ps %>% 
		mutate(p = round(p *100, round_digits)) %>%
		pivot_longer(-c(sex, age_cat_coarse, comm_type, type)) %>%
		mutate(
			sex = paste('sex_', sex, sep=''),
			age_cat_coarse = paste('age_cat_coarse_', age_cat_coarse, sep=''),
			comm_type = paste('comm_type_', comm_type, sep=''),
			type = case_when(
				type == 'participant-visits' ~ 'viremic_participant_visits',
				type == 'phsc' ~ 'phsc_participants')) %>%
		select(age_cat_coarse, sex, comm_type, name, type, value) %>%
		unite("label", age_cat_coarse:type)
	out[["labels"]] = c(labels, ps_dat$label)
	out[["vals"]] = c(vals, ps_dat$value)
	return(out)
}

calc_phsc_output = function(out, phsc_dat, round_digits=2){
	# unpack input
	labels = out[["labels"]]
	vals = out[['vals']]
	#### PHYLOSCANNER OUTPUT ####
	#labels = c(labels, "median_cov_all_windows")
	#vals = c(vals, median((phsc_dat %>% group_by(study_id, round) %>% summarise(n=n()))$n))
	#labels = c(labels, "sd_cov_all_windows")
	#vals = c(vals, round(sd((phsc_dat %>% group_by(study_id, round) %>% summarise(n=n()))$n), 2))
	all_n_d_i = summarise_n_d(tabulate_n_d(phsc_dat %>% filter(id_subgraph_reads > 0))) %>% 
		select(id, study_id, run, log10_copies, sequencing_technology, N_obs, MI_obs)
	uniq_n_d_i = summarise_n_d(tabulate_n_d(phsc_dat %>% filter(id_subgraph_reads > 0 & window_type == "unique"))) %>% 
		select(id, study_id, run, log10_copies, sequencing_technology, N_obs, MI_obs)
	n_d_i = all_n_d_i %>% full_join(uniq_n_d_i, by=c('id', 'study_id', 'run', 'log10_copies')) %>%
		mutate(N_obs.y = replace_na(N_obs.y, 0),
			MI_obs.y = replace_na(MI_obs.y, 0))
	labels=c(labels, 'all_unique_n_corr')
	vals=c(vals, round(cor(n_d_i$N_obs.x, n_d_i$N_obs.y, method = 'pearson'), round_digits))
	labels=c(labels, 'all_unique_mi_corr')
	vals=c(vals, round(cor(n_d_i$MI_obs.x, n_d_i$MI_obs.y, method = 'pearson'), round_digits))
	alt_uniq_n_d_i = summarise_n_d(tabulate_n_d(phsc_dat %>% filter(id_subgraph_reads > 0 & window_type == "unique_alt"))) %>% 
		select(id, study_id, run, log10_copies, sequencing_technology, N_obs, MI_obs)
	alt_n_d_i = uniq_n_d_i %>% full_join(alt_uniq_n_d_i, by=c('id', 'study_id', 'run', 'log10_copies')) %>%
		mutate(
			N_obs.x = replace_na(N_obs.x, 0),
			MI_obs.x = replace_na(MI_obs.x, 0),
			N_obs.y = replace_na(N_obs.y, 0),
			MI_obs.y = replace_na(MI_obs.y, 0))
	labels=c(labels, 'unique_unique_alt_n_corr')
	vals=c(vals, round(cor(alt_n_d_i$N_obs.x, alt_n_d_i$N_obs.y, method = 'pearson'), round_digits))
	labels=c(labels, 'unique_unique_alt_mi_corr')
	vals=c(vals, round(cor(alt_n_d_i$MI_obs.x, alt_n_d_i$MI_obs.y, method = 'pearson'), round_digits))
	labels = c(labels, 'n_unique_N_obs_gt_0')
	# observed windows with sequence data
	obs = uniq_n_d_i %>% filter(N_obs > 0)
	n_unique_n_obs_gt_0 = nrow(obs)
	vals = c(vals, format(n_unique_n_obs_gt_0, big.mark=','))
	labels = c(labels, 'p_unique_N_obs_gt_0')
	n_phsc_samples = as.numeric(str_replace(vals[which(labels == 'n_phsc_samples')],',',''))
	vals = c(vals, round(100*n_unique_n_obs_gt_0 / n_phsc_samples, round_digits))
	labels = c(labels, 'median_unique_N_obs')
	vals = c(vals, median(obs$N_obs))
	labels=c(labels, 'iqr_unique_N_obs')
	vals = c(vals, round(IQR(obs$N_obs), round_digits))
	window_size = phsc_dat$window_end[1] - phsc_dat$window_start[1] + 1
	labels = c(labels, 'median_unique_N_obs_bp')
	vals = c(vals, median(obs$N_obs*window_size))
	labels=c(labels, 'iqr_unique_N_obs_bp')
	vals = c(vals, round(IQR(obs$N_obs*window_size), round_digits))
	# stratified by viral load 
	labels = c(labels, 'median_unique_N_obs_3_3.5')
	vals = c(vals, median((obs %>% filter(log10_copies >= 3 & log10_copies < 3.5))$N_obs))
	labels = c(labels, 'iqr_unique_N_obs_3_3.5')
	vals = c(vals, round(IQR((obs %>% filter(log10_copies >= 3 & log10_copies < 3.5))$N_obs), round_digits))
	labels = c(labels, 'median_unique_N_obs_5_inf')
	vals = c(vals, median((obs %>% filter(log10_copies >= 5))$N_obs))
	labels = c(labels, 'iqr_unique_N_obs_5_inf')
	vals = c(vals, round(IQR((obs %>% filter(log10_copies >= 5))$N_obs), round_digits))
	n_unique_mi_obs_gt_0 = nrow(obs %>% filter(MI_obs > 0))
	# stratified y sequencing technology
	labels = c(labels, 'median_unique_N_obs_bc')
	vals = c(vals, median((obs %>% filter(sequencing_technology == 'bait_capture'))$N_obs))
	labels = c(labels, 'iqr_unique_N_obs_bc')
	vals = c(vals, round(IQR((obs %>% filter(sequencing_technology == 'bait_capture'))$N_obs), round_digits))
	labels = c(labels, 'median_unique_N_obs_amp')
	vals = c(vals, median((obs %>% filter(sequencing_technology == 'amplicon'))$N_obs))
	labels = c(labels, 'iqr_unique_N_obs_amp')
	vals = c(vals, round(IQR((obs %>% filter(sequencing_technology == 'amplicon'))$N_obs), round_digits))
	# how many with >0 multiple subgraph windows
	labels = c(labels, 'n_unique_MI_obs_gt_0')
	vals = c(vals, n_unique_mi_obs_gt_0)
	labels = c(labels, 'p_unique_MI_obs_gt_0')
	vals = c(vals, round(100*n_unique_mi_obs_gt_0 / n_unique_n_obs_gt_0, 2))
	labels = c(labels, 'p_unique_MI_obs_0')
	vals = c(vals, round(100*(n_unique_n_obs_gt_0 - n_unique_mi_obs_gt_0) / n_unique_n_obs_gt_0, 2))
	n_unique_mi_obs_gt_0_n_obs_29 = nrow(obs %>% filter(MI_obs > 0 & N_obs == 29))
	n_unique_n_obs_29 = nrow(obs %>% filter(N_obs == 29))
	labels = c(labels, 'n_unique_MI_obs_gt_0_N_obs_29')
	vals = c(vals, n_unique_mi_obs_gt_0_n_obs_29)
	labels = c(labels, 'n_unique_N_obs_29')
	vals = c(vals, nrow(obs %>% filter(N_obs == 29)))
	labels = c(labels, 'p_unique_MI_obs_gt_0_N_obs_29')
	vals = c(vals, round(100*n_unique_mi_obs_gt_0_n_obs_29 / n_unique_n_obs_29, round_digits))
	n_unique_mi_obs_gt_0_n_obs_lt_29 = nrow(obs %>% filter(MI_obs > 0 & N_obs < 29))
	n_unique_n_obs_lt_29 = nrow(obs %>% filter(N_obs < 29))
	labels = c(labels, 'n_unique_MI_obs_gt_0_N_obs_lt_29')
	vals = c(vals, n_unique_mi_obs_gt_0_n_obs_lt_29)
	labels = c(labels, 'n_unique_N_obs_lt_29')
	vals = c(vals, nrow(obs %>% filter(N_obs < 29)))
	labels = c(labels, 'p_unique_MI_obs_gt_0_N_obs_lt_29')
	vals = c(vals, round(100*n_unique_mi_obs_gt_0_n_obs_lt_29 / n_unique_n_obs_lt_29, round_digits))
	labels = c(labels, 'n_unique_participants_N_obs_gt_0')
	vals = c(vals, 
		format(nrow(uniq_n_d_i %>% 
			filter(N_obs > 0) %>% 
			select(study_id) %>% unique()), big.mark=','))
	labels = c(labels, 'n_unique_MI_obs_median')
	vals = c(vals, median((obs %>% filter(MI_obs > 0))$MI_obs))
	labels = c(labels, 'n_unique_MI_obs_iqr')
	vals = c(vals, IQR((obs %>% filter(MI_obs > 0))$MI_obs))
	labels = c(labels, 'p_unique_N_obs_w_MI_obs_median')
	vals = c(vals, round(
		median(100*(obs %>% filter(MI_obs > 0))$MI_obs/
				(obs %>% filter(MI_obs > 0))$N_obs),
		round_digits))
	labels = c(labels, 'p_unique_N_obs_w_MI_obs_iqr')
	vals = c(vals, round(
		IQR(100*(obs %>% filter(MI_obs > 0))$MI_obs/
				(obs %>% filter(MI_obs > 0))$N_obs),
		round_digits))
	# cophenetic distances
	x = phsc_dat %>% 
		filter(window_type == "unique") %>%
		select(id, window_start, id_min_sg_dist) %>%
		filter(!is.na(id_min_sg_dist)) %>%
		unique()
	labels = c(labels, 'median_multi_subgraph_dist')
	vals = c(vals, round(median(x$id_min_sg_dist), round_digits))
	labels = c(labels, 'iqr_multi_subgraph_dist')
	vals = c(vals, round(IQR(x$id_min_sg_dist), round_digits))
	# sexpever responses among men, including only first participant visits
	men_first_pv = summarise_n_d(tabulate_n_d(
			phsc_dat %>% filter(id_subgraph_reads > 0 & window_type == 'unique'))) %>% 
		select(study_id, sex, sexpever, plhiv_sexpever_std, round, N_obs, log10_copies) %>%
		group_by(study_id) %>%
		filter(N_obs == max(N_obs)) %>%
		filter(log10_copies == max(log10_copies)) %>%
		filter(sex == 'M')
	labels = c(labels, 'n_men_par_sexpeverNum')
	vals = c(vals, nrow(men_first_pv %>% 
		filter(plhiv_sexpever_std != 93 & !is.na(sexpever) & 
			plhiv_sexpever_std != 100)))
	labels = c(labels, 'p_men_par_sexpeverNum')
	vals = c(vals, 
		round(100*nrow(men_first_pv %>% 
				filter(plhiv_sexpever_std != 93 & !is.na(plhiv_sexpever_std) & 
					plhiv_sexpever_std != 100))/
			nrow(men_first_pv), round_digits))
	labels = c(labels, 'n_men_par_sexpeverMany')
	vals = c(vals, nrow(men_first_pv %>% filter(plhiv_sexpever_std == 93)))
	labels = c(labels, 'p_men_par_sexpeverMany')
	vals = c(vals, 
		round(100*nrow(men_first_pv %>% filter(plhiv_sexpever_std == 93))/
			nrow(men_first_pv), round_digits))
	labels = c(labels, 'n_men_par_sexpeverMissing')
	vals = c(vals, nrow(men_first_pv %>% filter(plhiv_sexpever_std == 100)))
	labels = c(labels, 'p_men_par_sexpeverMissing')
	vals = c(vals, 
		round(100*nrow(men_first_pv %>% filter(plhiv_sexpever_std == 100))/
			nrow(men_first_pv), round_digits))
	# repack output
	out[["labels"]] = labels
	out[["vals"]] = vals
	return(out)
}


calc_fit_output = function(out, name, param, f){
	# unpack input
	labels = out[["labels"]]
	vals = out[['vals']]
	if (!all(is.na(param))){
		# tabulate number of simulated samples
		labels = c(labels, paste(c('n_', name, '_simulated_samples'), collapse=''))
		vals = c(vals, format(as.numeric((param %>% filter(arg == 'N_ind'))$value), big.mark=','))
	}
	# rates so multiply by 100
	for (r in c(
			# average prob MI across all samples
			# fpr & fnr
			f$name[(grepl('prob_mi',f$name) & !grepl('\\[', f$name) & !grepl('coeffs', f$name))| 
				(f$name == 'prob_mi_fpr') | 
				(f$name == 'prob_mi_fnr')],
			# strata prob mi captured in the above
			# weighted prob MI from poststratification, includes strata weighted prob MI
			f$name[grepl('wprob_mi', f$name)])){
		r_out = str_replace(str_replace(r, '\\]', ''), '\\[', '')
		# simulated value, if included in parameters
		if (!all(is.na(param))){
			if (r %in% param$arg){
				# simulated FPR
				labels = c(labels, paste(c(name, '_', r_out, '_percent'), collapse=''))
				vals = c(vals, as.numeric((param %>% filter(arg == r))$value)*100)
			}
		}
		# fit value, if included in fit
		f_r = f %>% filter(name == r)
		for (col in c('median', 'lower', 'upper')){
			labels = c(labels, paste(c(name, '_fit_', r_out, '_', col, '_percent'), collapse=''))
			vals = c(vals, round(100*f_r[col][[1]], round_digits))
		}
		col = 'median'
		if (r == 'prob_mi' | grepl('wprob_mi', r)){
			labels = c(labels, paste(c(name, '_fit_', r_out, '_', col, '_percent_int'), collapse=''))
			vals = c(vals, round(100*f_r[col][[1]], 0))
		}
	}
	# not rates so don't multiply 
	for (r in c(
			f$name[grepl("logit_prob_seq", f$name) & !grepl("ind", f$name)],
			f$name[grepl('logit_prob_seq_ind_sd', f$name)], 
			f$name[grepl("logit_prob_mi_baseline", f$name)],
			f$name[grepl("logit_prob_mi_fpr", f$name)],
			f$name[grepl("logit_prob_mi_fnr", f$name)],
			f$name[grepl("logit_prob_mi_coeffs", f$name)],
		  	f$name[grepl('tau', f$name)], 
		  	f$name[grepl('lambda', f$name)],
		  	# strata prevance ratio
			f$name[grepl('pr_mi_', f$name)],
			# strata multivariate risk ratio
			f$name[grepl('multivar_rr_mi_', f$name)],
			# weighted strava prevalence ratio from poststratification
			f$name[grepl('wpr_mi_', f$name)])){
		r_out = str_replace(str_replace(r, '\\]', ''), '\\[', '')
		# if there are simulation parameters provided
		if (!all(is.na(param))){
			if (r %in% param$arg){
				# simulated FPR
				labels = c(labels, paste(c(name, '_', r_out), collapse=''))
				vals = c(vals, round(as.numeric((param %>% filter(arg == r))$value), round_digits))
			}
		}
		f_r = f %>% filter(name == r)
		r_round_digits = round_digits
		if (r == "prob_mi_fpr"){r_round_digits = round_digits*2}
		for (col in c('median', 'lower', 'upper')){
			labels = c(labels, paste(c(name, '_fit_', r_out, '_', col), collapse=''))
			vals = c(vals, round(f_r[col][[1]], r_round_digits))
		}
		for (col in c('bulk_ess', 'tail_ess', 'rhat')){
			labels = c(labels, paste(c(name, '_fit_', r_out, '_', col), collapse=''))
			vals = c(vals, round(f_r[col][[1]], round_digits))
		}
	}
	if (any(!is.na(f$idx))){
		labels = c(labels, paste(name, '_n_participant_visits', sep=''))
		vals = c(vals, format(f %>% select(idx) %>% drop_na() %>% max(), big.mark=','))
	}
	# repack output
	out[["labels"]] = labels
	out[["vals"]] = vals
	return(out)
}


#### MAIN FUNCTION ####
p <- arg_parser("plot summary of simulated data")
p <- add_argument(p, "--phscDat", help="phyloscanner data file", nargs=1)
p <- add_argument(p, "--metadata", help="rccs data file", nargs=1)
p <- add_argument(p, "--fits", help="model fit summaries", nargs=Inf)
p <- add_argument(p, "--labels", help="model fit labels", nargs=Inf)
p <- add_argument(p, "--params", help="model fit params", nargs=Inf)
p <- add_argument(p, "--empiricalCutoffFit", help="model fit params", nargs=Inf)
p <- add_argument(p, "--empiricalSexpeverMenFit", help="model fit params", nargs=Inf)
args <- parse_args(p)

round_digits = 2

#args$phscDat = 'output/211220_allreads_phsc_all_subgraphs_format_par.tsv' 
#args$metadata = 'output/211220_allreads_phsc_metadata.tsv'
#args$fits =c(
#	'fit/base_simulation_base_deep-phyloMI__summary.tsv',
#	'fit/full_simulation_base_deep-phyloMI__summary.tsv',
#	'fit/full_simulation_deep-phyloMI__summary.tsv',
#	'fit/extended_simulation_deep-phyloMI__summary.tsv',
#	'fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_age_sex_comm_summary.tsv',
#	'fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_age_sex_comm_alt_summary.tsv',
#	'fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_age_sex_comm_incl_minor_summary.tsv',
#	'fit/211220_allreads_phsc_all_subgraphs_format_par_amplicon_tmp_age_sex_comm_amplicon_summary.tsv',
#	'fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_seq_summary.tsv',
#	'fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_comm_seq_summary.tsv',
#	'fit/211220_allreads_phsc_all_subgraphs_format_par_m_deep-phyloMI_sexpever_men_summary.tsv', 
#	'fit/211220_allreads_phsc_all_subgraphs_format_par_m_complete_deep-phyloMI_sexpever_men_complete_summary.tsv',
#	'fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_var_select_summary.tsv')

#args$labels = c(
#	'base_base',
#	'full_base',
#	'full_full',
#	'ext_ext',
#	'empirical_age_sex_comm',
#	'empirical_age_sex_comm_alt',
#	'empirical_age_sex_comm_incl_minor',
#	'empirical_age_sex_comm_amplicon',
#	'empirical_seq',
#	'empirical_comm_seq',
#	'empirical_sexpever_men',
#	'empirical_sexpever_men_complete',
#	'empirical_var_select')
#args$params = c(
#	'simulations/base_simulation_params.tsv',
#	'simulations/full_simulation_params.tsv',
#	'simulations/full_simulation_params.tsv',
#	'simulations/extended_simulation_params.tsv',
#	'NA',
#	'NA',
#	'NA',
#	'NA',
#	'NA',
#	'NA',
#	'NA')
#args$empiricalCutoffFit = 'fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_age_sex_comm_summary.tsv'
#args$empiricalSexpeverMenFit = 'fit/211220_allreads_phsc_all_subgraphs_format_par_m_deep-phyloMI_sexpever_men_summary.tsv'


args$params = if_else(args$params == 'NA', NA, args$params)

# read in data
phsc_dat = read_tsv(args$phscDat, show_col_types=FALSE)
metadata = read_tsv(args$metadata, show_col_types=FALSE) 

out = list()
out[["labels"]] = c()
out[["vals"]] = c()

# dataset statistics
out = calc_basic_stats(out, metadata)
out = calc_epi_stats(out, metadata)
out = calc_supp_epi_stats(out, metadata, round_digits=2)
out = calc_vl_stats(out, metadata, round_digits=round_digits)
out = calc_phsc_input(out, metadata)
out = calc_phsc_config(out, phsc_dat)
out = calc_weight_stats(out, metadata, round_digits)

# %>%
#		add_row(tibble_row(
#			arg='logit_prob_mi_fpr', 
#			value=as.character(round(logit(as.numeric(.$value[.$arg == "prob_mi_fpr"])), round_digits)))) %>% 
#		add_row(tibble_row(
#			arg='logit_prob_mi_fnr', 
#			value=as.character(round(logit(as.numeric(.$value[.$arg == "prob_mi_fnr"])), round_digits))))

# model parameters and results
for (fit_idx in 1:length(args$fits)){
	if (!is.na(args$params[fit_idx])){
		params = read_tsv(args$params[fit_idx], show_col_types=FALSE) %>%
		add_row(tibble_row(
			arg='logit_prob_mi_fpr', 
			value=as.character(round(logit(as.numeric(.$value[.$arg == "prob_mi_fpr"])), round_digits)))) %>% 
		add_row(tibble_row(
			arg='logit_prob_mi_fnr', 
			value=as.character(round(logit(as.numeric(.$value[.$arg == "prob_mi_fnr"])), round_digits))))
	}else{
		params = NA
	}
	out = calc_fit_output(out, args$labels[fit_idx], params, 
		read_tsv(args$fits[fit_idx], show_col_types=FALSE))
	# check for post-stratification
	if (file.exists(gsub("_summary.tsv", "_poststrat_summary.tsv", args$fits[fit_idx]))){
		out = calc_fit_output(out, args$labels[fit_idx], params, 
		read_tsv(gsub("_summary.tsv", "_poststrat_summary.tsv", args$fits[fit_idx]), show_col_types=FALSE) %>% mutate(idx=NA))
	}
}

out = calc_phsc_output(out, phsc_dat, round_digits=round_digits)

# unpack
vals = out[["vals"]]
labels = out[["labels"]]

# add amplicon stats
amp_dat = summarise_n_d(
	read_tsv(gsub(".tsv", "_amplicon.tsv", args$phscDat), show_col_types=FALSE))
n_unique_n_obs_gt_0_amplicon = nrow(amp_dat %>% filter(N_obs > 0))
n_unique_mi_obs_gt_0_amplicon = nrow(amp_dat %>% filter(MI_obs > 0))
labels = c(labels, 'n_unique_N_obs_gt_0_amplicon')
vals = c(vals, format(n_unique_n_obs_gt_0_amplicon, big.mark=','))
labels = c(labels, 'n_unique_MI_obs_gt_0_amplicon')
vals = c(vals, n_unique_mi_obs_gt_0_amplicon)
labels = c(labels, 'p_unique_MI_obs_gt_0_amplicon')
vals = c(vals, round(100*n_unique_mi_obs_gt_0_amplicon / n_unique_n_obs_gt_0_amplicon, 2))


# additional empirical results
# proportion of sequenced participants by sequencing technology and community type
par_dat = read_tsv(args$phscDat, show_col_types=FALSE)
p_fishing = par_dat %>% select(study_id, round, comm_type, sequencing_technology) %>%
	unique() %>% group_by(sequencing_technology) %>% 
	summarise(
		n=n(),
		k_fishing = sum(comm_type == "fishing"),
		p_fishing = k_fishing/n)

labels = c(labels, 'p_amplicon_fishing')
vals = c(vals, round((p_fishing %>% filter(sequencing_technology == "amplicon"))$p_fishing*100, round_digits))

labels = c(labels, 'p_bait_fishing')
vals = c(vals, round((p_fishing %>% filter(sequencing_technology == "bait_capture"))$p_fishing*100, round_digits))

# compare multiple subgraph mrca dists to all mrca dists

#all_mrca_dists = phsc_dat %>% select(study_id, round, window_start, mrca_dist) %>% unique()
#median_mi_mrca_dist = median((phsc_dat %>% select(study_id, round, window_start, dist) %>% filter(!is.na(dist)) %>% unique())$dist)
#labels = c(labels, 'all_mrca_dist_median_mi_dist_percentile'))
#vals = c(vals, round(100*sum(all_mrca_dists$mrca_dist >= median_mi_mrca_dist)/nrow(all_mrca_dists), 
#	round_digits)

# posterior n_mi with prob_mi > 50
#args$empiricalCutoffFit = 'fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_age_sex_comm_summary.tsv'
print(args$empiricalCutoffFit)
n_mi_cutoff_dat = read_tsv(str_replace(args$empiricalCutoffFit, '_summary.tsv', '_n_mi_dichotomized.tsv'), 
	show_col_types=FALSE)


labels = c(labels, 'empirical_age_sex_comm_fit_ind_prob_mi_ge_50_median')
vals = c(vals, (n_mi_cutoff_dat %>% filter(cutoff == "0.5"))$`0.50`)
labels = c(labels, 'empirical_age_sex_comm_fit_ind_prob_mi_ge_50_median_percent')
vals = c(vals, 
	round(100*((n_mi_cutoff_dat %>% filter(cutoff == "0.5"))$`0.50` / 
			length(unique((read_tsv(args$empiricalCutoffFit, show_col_types=FALSE) %>% filter(!is.na(study_id)))$study_id))),
		round_digits))


# MI participants to match posterior pred
labels = c(labels, 'empirical_age_sex_comm_fit_ind_prob_mi_ge_empirical_median')
vals = c(vals, format((n_mi_cutoff_dat %>% filter(cutoff == "n_mi"))$`0.50`, big.mark=','))
labels = c(labels, 'empirical_age_sex_comm_fit_ind_prob_mi_empirical_median_percent')
vals = c(vals, round(100*(n_mi_cutoff_dat %>% filter(cutoff == "empirical"))$`0.50`, round_digits))

# compare risk of MI based on sexpever among those in fishing communities
# need to get standardization value!
sexpever_meanFishing24_29 = (metadata %>% filter(finalhiv == 'P' & sex == 'M' & comm_type == 'fishing' &
	age_cat_fine == '(24,29]') %>%
	select(plhiv_sexpever_mean) %>% unique())[[1]]

sexpever_fit = readRDS(str_replace(args$empiricalSexpeverMenFit, "_summary.tsv", ".Rds"))
sexpever_draws = as_tibble(sexpever_fit$draws(
		variables = c('logit_prob_mi_baseline', 'logit_prob_mi_coeffs[1]', 'logit_prob_mi_coeffs[3]', 'logit_prob_mi_coeffs[4]'),
        inc_warmup = FALSE,
        format = "draws_df")) %>%
	select(-.chain, -.iteration, -.draw) %>%
	mutate(cit = seq(1,n()))


# risk of MI with 1 lifetime sexpartner
# THIS IS WRONG
rr_sexpever_Fishing24_29_1v30 = sexpever_draws %>%
	mutate(rr = inv_logit(
			logit_prob_mi_baseline + 
				`logit_prob_mi_coeffs[1]` + 
				`logit_prob_mi_coeffs[3]` * (30 - sexpever_meanFishing24_29) +
				`logit_prob_mi_coeffs[4]` * (30 - sexpever_meanFishing24_29)) / 
			inv_logit(
			logit_prob_mi_baseline + 
				`logit_prob_mi_coeffs[1]` + 
				`logit_prob_mi_coeffs[3]` * (1 - sexpever_meanFishing24_29) + 
				`logit_prob_mi_coeffs[4]` * (1 - sexpever_meanFishing24_29))) %>%
	summarise(median=median(rr), lower = get_hpd(rr)$lower, upper=get_hpd(rr)$upper)


for (col in colnames(rr_sexpever_Fishing24_29_1v30)){
	labels = c(labels, paste(c("empirical_sexpever_men_fit_rr_fishing2429_30v1_", col), collapse=''))
	vals = c(vals, round(rr_sexpever_Fishing24_29_1v30[[col]][1], round_digits))
}

out = tibble(labels=labels, vals=vals) %>% mutate(labels = gsub("&", "\\\\\\\\&", labels))

write_csv(out, 'output/statistics.csv', quote = "all")
# need this to manually replace values in tex file for pandoc
write_delim(out, 'output/statistics.scsv', delim=";")

