suppressMessages(require(tidyverse))
suppressMessages(require(haven))
suppressMessages(require(stringr))
suppressMessages(require(fitdistrplus))
suppressMessages(require(conflicted))
suppressMessages(source('scripts/utils.R'))
suppressMessages(conflict_prefer("select", "dplyr"))
suppressMessages(conflict_prefer("filter", "dplyr"))


#### ------------------------- ####
#### 0. DESIRED OUTPUT COLUMNS ####
#### ------------------------- ####
rccs_vars = c(
	'study_id', 'round', 'visit_dt',
	'finalhiv_orig', 'finalhiv', 'log10_copies', 'log10_copies_cat',
	'sex', 'age_cat_coarse', 'age_cat_fine', 
	'comm_type', 
	'male_circumcision', 
	'married', 
	'sexpever', 'plhiv_sexpeverSome_mean', 'plhiv_sexpeverMany_shape', 'plhiv_sexpeverMany_scale',
	'plhiv_sexpeverMissing_shape', 'plhiv_sexpeverMissing_scale',
	'plhiv_sexpever_impute_naive', 'plhiv_sexpever_mean', 'plhiv_sexpever_std',
	'barworker', 'in_migrant')

pg_vars = c("pangea_id", "id", "study_id", "visit_dt",  'sequencing_technology')


#### ----------------------------- ####
#### 1. FORMAT PHSC INPUT METADATA ####
#### ----------------------------- ####
# maps rename_id (AID*) to pangea_id
input = read_tsv('data/210120_RCCSUVRI_phscinput_samples.tsv', show_col_types=FALSE) %>%
	rename(
		id = RENAME_ID, 
		sample_id = SAMPLE_ID,
		unit_id = UNIT_ID) %>%
	mutate(
		study_id = str_split(unit_id, '-', simplify=TRUE)[,2],
		cohort = str_split(unit_id, "-", simplify=T)[,1],
		pangea_id = str_split(PANGEA_ID, "_", simplify=T)[,2],
		sequencing_technology = case_when(
			str_detect(sample_id, 'miseq') ~ 'amplicon',
			str_detect(sample_id, 'hiseq') ~ 'amplicon',
			!str_detect(sample_id, 'miseq') & !str_detect(sample_id, 'hiseq') ~ 'bait_capture'))


#### ------------------------- ####
#### 2. FORMAT PANGEA METADATA ####
#### ------------------------- ####
# this is only needed to map pangea_id to visit_dt
# do NOT pull copies from here as some don't match the visit_dt
pg_dat = read_csv('data/200316_pangea_db_sharing_extract_rakai.csv', show_col_types=FALSE)

# confirm only one pangea_id per participant visit
stopifnot(nrow(pg_dat %>% group_by(pangea_id) %>% 
	summarise(n=length(unique(visit_dt))) %>% filter(n == 2)) == 0)

# merge and get just rakai samples
d = input %>% left_join(pg_dat, by='pangea_id') %>% 
	filter(cohort == 'RK') %>% 
	# when resequenced, get one with highest readnumber
	group_by(study_id, visit_dt) %>%
	filter(readnum_hiv == max(readnum_hiv)) %>%
	select(all_of(pg_vars)) %>% unique()


#### ----------------------- ####
#### 3. FORMAT RCCS METADATA ####
#### ----------------------- ####
# read in rccs data
rccs_dat = read_format_rccs_dat('data/RCCSdata_R001_R019_VOIs.dta')
# filter for inclustion criteria
rccs_dat = rccs_dat %>% 
	filter(ageyrs >= 15 & ageyrs <= 49 & round != 15.1)
# read in full occupation data
occup = read_csv('data/quest_occup1_occup2_R1419.csv', show_col_types = FALSE) %>%
	mutate(round = as.numeric(str_replace(str_remove(round, "^R0+"), 'S', '.1')))
barworker = apply(occup[,c("occup1", "occup2")] == 12 | 
	occup[,c("occup1", "occup2")] == 23 | occup[,c("occup1", "occup2")] == 18, 1, any)
check_missing = function(x){
	return(is.na(x) | x == 97 | x == 98)
}
barworker = if_else(apply(check_missing(occup[,c("occup1", "occup2")]), 1, all), NA, barworker)
barworker = tibble(study_id = occup$study_id, round=occup$round, barworker_occup = barworker)

# add in partner occupation data 
partner_occup = read_csv('data/partnerblock_occup_R1419.csv', show_col_types = FALSE) %>%
	mutate(round = as.numeric(str_replace(str_remove(round, "^R0+"), 'S', '.1')))
# for each participant, do they report having sex with any barworkers
# a bit naive as not done on per-partner level, todo, fix
sex_w_barworker = apply(partner_occup[,3:ncol(partner_occup)] == 12 | 
	partner_occup[,3:ncol(partner_occup)] == 23 | partner_occup[,3:ncol(partner_occup)] == 18, 1, any)
check_missing = function(x){
	return(is.na(x) | x == 97 | x == 98)
}
sex_w_barworker = if_else(apply(check_missing(partner_occup[,3:ncol(partner_occup)]), 1, all), NA, sex_w_barworker)
sex_w_barworker = tibble(study_id = partner_occup$study_id, round=partner_occup$round, sex_w_barworker = sex_w_barworker)

# sexyear == 2 -> no sex in last year,therefore no sex w/ barworker in past year
# sexp1yr == 0 -> no sex in last year,therefore no sex w/ barworker in past year
rccs_dat = rccs_dat %>% 
	left_join(barworker, 
		by=c('study_id', 'round')) %>%
	left_join(sex_w_barworker, 
		by=c('study_id', 'round')) %>%
	mutate(sex_w_barworker_final = case_when(
		sex_w_barworker == TRUE ~ TRUE,
		sex_w_barworker == FALSE ~ FALSE,
		sexyear == 2 ~ FALSE, 
		sexp1yr == 0 ~ FALSE))
	
# format columns
r = rccs_dat %>%
	mutate(
		finalhiv_orig = if_else(finalhiv == "", NA, finalhiv),
		finalhiv = if_else(finalhiv == "P", "P", "N"),
		visit_dt = int_date, 
		age_cat_coarse = cut(ageyrs, breaks=c(14, 24, 34, 49)),
		age_cat_fine = cut(ageyrs, breaks=seq(14,49,5)),
		married = case_when(
			currmarr == 2 ~ FALSE,
			currmarr == 1 ~ TRUE,
			evermarr == 2 ~ FALSE),
		#barworker_old = case_when(
		#	sex == 'F' & occup1 == 12 ~ TRUE,
		#	sex == 'F' & occup1 == 23 ~ TRUE,
		#	sex == 'F' & occup1 != 12 & occup1 != 23 & occup1 != 97 & occup1 != 98 & !is.na(occup1) ~ FALSE, 
		#	sex == 'M' & sex_w_barworker_final == TRUE ~ TRUE,
		#	sex == 'M' & sex_w_barworker_final == FALSE ~ FALSE),
		barworker = case_when(
			sex == 'F' & barworker_occup == TRUE ~ TRUE,
			sex == 'F' & barworker_occup == FALSE ~ FALSE,
			sex == 'F' & is.na(barworker_occup) & !is.na(occup1) & (occup1 == 12 | occup1 == 23) ~ TRUE,
			sex == 'F' & is.na(barworker_occup) & !is.na(occup1) & (occup1 != 12 & occup1 != 23) ~ FALSE,
			sex == 'M' & sex_w_barworker_final == TRUE ~ TRUE,
			sex == 'M' & sex_w_barworker_final == FALSE ~ FALSE),
		comm_type = if_else(
			as_factor(comm_type) == 'Fishing', 'fishing', 'inland'),
		in_migrant = case_when(
			mobility == 3 ~ TRUE,
			mobility == 9 ~ NA,
			is.na(mobility) ~ NA,
			mobility != 3 & mobility != 9 & !is.na(mobility) ~ FALSE),
		male_circumcision = case_when(
			sex == 'M' & circum == 1 ~ TRUE,
			sex == 'M' & circum == 2 ~ FALSE,
			sex == 'M' & !(circum == 1 | circum == 2) ~ NA,
			sex != 'M' ~ FALSE),
		#sexpever = if_else(finalhiv == 'P' & sexpever == 0, NA, sexpever),
		#arv = cuarvmed == '1' | arvmed == '1'
		)	


# MEAN SEXPEVER VALUES WITHIN STRATA
group_vars = c('sex', 'age_cat_fine', 'comm_type')
# at present, the categorical responses we consider are:
# sexpever == 92: 1-2
# sexpever == 93: 3+
# sexpever == 0: missing
# sexpever == 97: missing/skip
# sexpever == 98: missing/skip
# sexpever == 99: missing
# these are hard coded, to do fix

# first, fit distribution to responses of 1+
plhiv_sexpeverMissing_fit = 
	bind_rows(r %>% filter(finalhiv == 'P' & 
			(!is.na(sexpever) & 
				sexpever >= 1 & 
				sexpever != 93 & 
				sexpever != 97 & 
				sexpever != 98 & 
				sexpever != 99)) %>%
		mutate(
			sexpever = as.numeric(sexpever),
			age_cat_fine = as.character(age_cat_fine)) %>%
		group_by_at(group_vars) %>%
		group_map(~c(unlist(.y), setNames(fitdist(.x$sexpever, 'lnorm')$estimate,
			c('plhiv_sexpeverMissing_shape', 'plhiv_sexpeverMissing_scale'))))) %>%
		mutate(plhiv_sexpeverMissing_shape = as.numeric(plhiv_sexpeverMissing_shape),
			plhiv_sexpeverMissing_scale = as.numeric(plhiv_sexpeverMissing_scale))

# next, mean value between response of 1 and 2
plhiv_sexpeverSome_mean = 
	r %>% filter(finalhiv == 'P' & (sexpever == 1 | sexpever == 2)) %>%
	group_by_at(group_vars) %>%
	summarise(plhiv_sexpeverSome_mean = mean(sexpever), .groups='drop')

# finally, fit distribution to responses of 3+
plhiv_sexpeverMany_fit = 
	bind_rows(r %>% filter(finalhiv == 'P' & 
			(!is.na(sexpever) & 
				sexpever >= 3 & 
				sexpever != 93 & 
				sexpever != 97 & 
				sexpever != 98 & 
				sexpever != 99)) %>%
		mutate(
			sexpever = as.numeric(sexpever),
			age_cat_fine = as.character(age_cat_fine)) %>%
		group_by_at(group_vars) %>%
		group_map(~c(unlist(.y), setNames(fitdist(.x$sexpever, 'lnorm')$estimate,
			c('plhiv_sexpeverMany_shape', 'plhiv_sexpeverMany_scale'))))) %>%
		mutate(plhiv_sexpeverMany_shape = as.numeric(plhiv_sexpeverMany_shape),
			plhiv_sexpeverMany_scale = as.numeric(plhiv_sexpeverMany_scale))


# merge in
# naive imputation, pulls means value for uncertain responses
# add mean sexpever
# assumes lognormal distribution!
r = r %>% 
	left_join(plhiv_sexpeverSome_mean, by=group_vars) %>%
	left_join(plhiv_sexpeverMany_fit, by=group_vars) %>%
	left_join(plhiv_sexpeverMissing_fit, group_vars) %>%
	mutate(
		plhiv_sexpever_impute_naive = case_when(
			finalhiv == 'P' & 
					(sexpever == 0  | 
						sexpever == 97 |
						sexpever == 98 |
						sexpever == 99) 
				~ exp(plhiv_sexpeverMissing_shape + 0.5*plhiv_sexpeverMissing_scale^2),
			finalhiv == 'P' & sexpever == 92 ~ plhiv_sexpeverSome_mean,
			finalhiv == 'P' & sexpever == 93 ~ exp(plhiv_sexpeverMany_shape + 0.5*plhiv_sexpeverMany_scale^2),
			finalhiv == 'P' & 
					sexpever != 0 & 
					sexpever != 97 & 
					sexpever != 98 & 
					sexpever != 99 & 
					sexpever != 92 & 
					sexpever != 93
				~ sexpever)) %>%
	group_by_at(group_vars) %>%
	mutate(
		plhiv_sexpever_mean = case_when(
			finalhiv == 'P' ~ mean(plhiv_sexpever_impute_naive, na.rm=TRUE)),
		plhiv_sexpever_std = case_when(
			finalhiv == 'P' & 
					(sexpever == 0  | 
						sexpever == 97 |
						sexpever == 98 |
						sexpever == 99)
				~ 100,								# taking as missing values
			finalhiv == 'P' & sexpever == 92 ~ 92,	# ambiguous response, 1-2
			finalhiv == 'P' & sexpever == 93 ~ 93,	# ambiguous response, 3+
			finalhiv == 'P' & 
					sexpever != 0 & 
					sexpever != 97 & 
					sexpever != 98 & 
					sexpever != 99 & 
					sexpever != 92 & 
					sexpever != 93
				~ sexpever - plhiv_sexpever_mean # standardize to mean value
			)) %>%
	ungroup()


#### ----------------------------------- ####
#### 4. MERGE INTO CONSOLIDATED METADATA ####
#### ----------------------------------- ####
# get the rounds that are present in the pangea data
get_rounds = unique((d %>% select(study_id, visit_dt) %>%
	left_join(r %>% select(study_id, round, visit_dt) , 
		by=c('study_id', 'visit_dt')) %>% drop_na())$round)

# merge into a single metadata file
metadata = r %>% 
	filter(round %in% get_rounds) %>%
	#mutate(n_comm_participant = length(unique(comm_num))) %>%
	select(all_of(rccs_vars)) %>%
	left_join(d %>% select(all_of(pg_vars)), 
		by=c('study_id', 'visit_dt'))	

# create new study_id
# tabulate number of communities
metadata = metadata %>% mutate(rccs_study_id = study_id) %>%
	group_by(study_id) %>%
	mutate(study_id = as.character(cur_group_id())) %>%
	ungroup() %>%
	mutate(
		study_id = str_pad(study_id, max(nchar(study_id)), pad='0', side='left'))

# get sampling dates
# sampling dates 
survey_dates = 
	metadata %>%
	# filter outlier dates
	group_by(round) %>%
	filter(visit_dt > quantile(visit_dt, c(5E-4),type=1) & visit_dt < quantile(visit_dt, c(1-5E-4),type=1)) %>%
	ungroup() %>%
	summarise(
		min = min(visit_dt, na.rm=TRUE), 
		max = max(visit_dt, na.rm=TRUE))

metadata$min_survey_date = paste(c(
				format(survey_dates$min, "%B"), 
				" ",
				format(survey_dates$min, "%Y")), collapse='')

metadata$max_survey_date = paste(c(
				format(survey_dates$max, "%B"), 
				" ",
				format(survey_dates$max, "%Y")), collapse='')

# add median round years
metadata = metadata %>% group_by(round) %>% 
	mutate(round_median_year = str_split(median(visit_dt), "-", simplify=T)[,1])


write_tsv(metadata, 'data/input_metadata_internal.tsv')
#write_tsv(metadata %>% select(-rccs_study_id, -comm_num, -visit_dt), 'data/input_metadata.tsv')




