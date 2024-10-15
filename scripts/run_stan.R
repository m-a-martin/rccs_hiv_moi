suppressMessages(require(cmdstanr))
suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(require(HDInterval))
suppressMessages(require(posterior))
suppressMessages(source('scripts/utils.R'))


scale_vars = function(scale_vars, d){
  # scales numeric variables to mean 0 and sd 1
  if (!all(is.na(scale_vars))){
    for (i in scale_vars){
      i_split = str_split(i, ":", simplify=TRUE)
      to_scale = i_split[1]
      if (length(i_split) > 1){
        group_vars_scale = str_split(i_split[2], ",", simplify=TRUE)[1,]
        d = d %>% group_by_at(group_vars_scale) %>%
          mutate(!!to_scale := scale(!!as.name(to_scale))[,1])
        }else{
          d = d %>% 
            mutate(!!to_scale := scale(!!as.name(to_scale))[,1])
        }
      print(d %>% summarise(
        m = mean(!!as.name(to_scale)), 
        s = sd(!!as.name(to_scale))) %>%
        mutate(var = to_scale))
      d = d %>% ungroup()
    }
  }
  return(d)
}


run_stan = function(args){
  #### OUTPUT OBJECTS ####
  out = list()
  out$mi_design_matrix_cols = c()
  out$seq_design_matrix_cols = c()
  #### READ IN AND PREPARE DATA ####
  d = read_tsv(args$dat, show_col_types=FALSE) 
  if ('round' %in% colnames(d)){
    d = d %>% mutate(round = as.character(round))
  }
  if (args$inputType == 'full'){
    tmp = d %>% 
      filter(window_type == args$getWindowType)
    d = summarise_n_d(tabulate_n_d(d %>% 
      filter(window_type == args$getWindowType)))
  }
  d = scale_vars(args$scaleVars, d)
  #### PREPARE STAN DATA ####
  adapt = 0.8
  stan_data = list()
  # first get all vars involved in both design matrices
  # and drop na values
  vars = c(
    unlist(str_split(unlist(str_split(
      args$seqDesignMatrix, ":", simplify=TRUE)), "\\*", simplify=TRUE)),
    unlist(str_split(unlist(str_split(
      args$miDesignMatrix, ":", simplify=TRUE)), "\\*", simplify=TRUE)))
  vars = vars[!is.na(vars) & vars != ""]
  d_f = d[apply(!is.na(d[,vars]),1,all),]
  n_rm = nrow(d)-nrow(df)
  if (length(n_rm) == 0){n_rm="0"}else{n_rm=as.character(n_rm)}
  print(paste(n_rm, 
    ' rows dropped due to NA values', sep=''))
  # seq design matrix
  if (!all(is.na(args$seqDesignMatrix))){
    # if data on sequencing technology, 
    # stratify sequencing success by technology
    seq_f =  paste('~',
      paste(args$seqDesignMatrix, collapse='+'), sep='')
    seq_design_matrix = as.matrix(model.matrix(
      as.formula(seq_f), data=d_f))
    out$seq_design_matrix_cols = colnames(seq_design_matrix)[2:ncol(seq_design_matrix)]
    stan_data$K_seq = as.integer(ncol(seq_design_matrix)-1)
    stan_data$X_seq = as.matrix(seq_design_matrix[,2:ncol(seq_design_matrix)])
  }else{
    # I'm not actually sure if this works
    # if there are no predictors
    # need to double check
    seq_design_matrix = vector('numeric')
    out$seq_design_matrix_cols = NA
    stan_data$K_seq = 0
    stan_data$X_seq = seq_design_matrix
  }
  # mi design matrix
  if (!all(is.na(args$miDesignMatrix))){
    # reorder factors for comm_type
    if ('comm_type' %in% colnames(d)){
      d_f$comm_type = factor(d_f$comm_type, levels = sort(unique(d_f$comm_type), decreasing=TRUE))
    }
    stan_data$nu = 2
    f = paste('~', paste(args$miDesignMatrix, collapse='+'), sep='')
    mi_design_matrix = as.matrix(model.matrix(
      as.formula(f), data=d_f))
    # standardize continous variables within strata 
    out$mi_design_matrix_cols = colnames(mi_design_matrix)[2:ncol(mi_design_matrix)]
    stan_data$K_mi_risk_factors = ncol(mi_design_matrix)-1
    stan_data$X_mi_risk_factors = as.matrix(mi_design_matrix[,2:ncol(mi_design_matrix)], 
      ncol=ncol(mi_design_matrix)-1)
    adapt = 0.95
  }
  # missing data 
  if (!all(is.na(args$miMissingDat))){
    idx_missing = c()
    missing_type = c()
    missing_prior = c()
    missing_std = c()
    missing_min = c()
    missing_max = c()
    for (idx in 1:length(args$miMissingDat)){
      i = args$miMissingDat[idx]
      i_split = str_split(i, ":", simplify=TRUE)
      attr_split = str_split(i_split[,2], ";", simplify=TRUE)
      missing_dat_col = i_split[,1]
      # WHICH VALUES NEED IMPUTING
      missing_dat_val = attr_split[,1]
      missing_data_min = attr_split[,2]
      missing_data_max = attr_split[,3]
      missing_data_shape = attr_split[,4]
      missing_data_scale = attr_split[,5]
      missing_data_std = attr_split[,6]
      where_missing = which(sweep(
        mi_design_matrix[,2:ncol(mi_design_matrix)] == missing_dat_val,
          MARGIN=2,
          grepl(missing_dat_col, colnames(mi_design_matrix)[2:ncol(mi_design_matrix)]),
          '*') == 1, 
        arr.ind=TRUE)
      idx_missing = rbind(idx_missing, where_missing)
      # MINIMUM FOR IMPUTED VALUES
      missing_min = c(missing_min, rep(as.numeric(missing_data_min), nrow(where_missing)))
      # MAXIMUM FOR IMPUTED VALUES 
      missing_max = c(missing_max, rep(as.numeric(missing_data_max), nrow(where_missing)))
      # PRIOR FOR IMPUTED VALUES
      missing_prior = rbind(missing_prior, 
        as.matrix(d_f[where_missing[,1],c(missing_data_shape, missing_data_scale)]))
      # STANDARDIZATION FOR IMPUTED VALUES
      if (length(attr_split)){
        if (tolower(missing_data_std) != 'none' & !is.na(missing_data_std)){
          missing_std = c(missing_std, d_f[where_missing[,1],][[missing_data_std]])
        }else{
          missing_std = c(missing_std, rep(0, nrow(where_missing)))
        }
      }else{missing_std = c(missing_std, rep(0, nrow(where_missing)))}
    }
    stan_data$N_missing = nrow(idx_missing)
    stan_data$idx_missing = idx_missing
    stan_data$X_mi_missing_std = missing_std
    stan_data$missing_prior = missing_prior
    # shared across all variables
    stan_data$missing_max = missing_max
    stan_data$missing_min = missing_min
 
    # format missing values appropriately 
  }else{
    # no missing values
    stan_data$N_missing = 0
    stan_data$idx_missing = matrix(0, 0,2)
    stan_data$X_mi_missing_std = vector('numeric')
    stan_data$missing_prior = vector('numeric')
    stan_data$missing_min = vector('numeric')
    stan_data$missing_max = vector('numeric')
  }
  stan_data$N_ind = nrow(d_f)
  stan_data$N_obs_max = d_f$N_windows[1]
  stan_data$N_obs = d_f$N_obs
  stan_data$MI_obs = d_f$MI_obs
  print("here")
  print(head(seq_design_matrix))
  print(head(stan_data$X_seq))
  print(stan_data$K_seq)
  #### COMPILE STAN MODEL ####
  model_compiled = cmdstanr::cmdstan_model(args$stan)
  #### RUN STAN ####
  out$fit = model_compiled$sample(
    data = stan_data,
    seed = 42,
    chains = 4,
    parallel_chains = 4,
    iter_warmup = 500,
    iter_sampling = 2000,
    refresh = 500, # print update every 500 iters,
    save_warmup = TRUE,
    show_messages = TRUE,
    show_exceptions = TRUE,
    adapt_delta = adapt)
  out$d_f = d_f
  out$stan_data = stan_data
  return(out)
}


#### --------------------- ####
#### MAIN CODE BEGINS HERE ####
#### --------------------- ####
p = arg_parser("run stan on input data")
#### INPUT DATA ####
p = add_argument(p, "--inputType", help="format of input file", default='summarized', type="character", nargs=1)
p = add_argument(p, "--getWindowType", help="which set of windows to get", default='unique', type="character", nargs=1)
p = add_argument(p, "--dat", help="input data file", nargs=1)
p = add_argument(p, "--stan", help="stan model file", nargs=1)
#### DATA TO SCALE ####
p = add_argument(p, "--scaleVars", 
  help="list of variables to scale. expected format is var:group_vars where group_vars is a semicolon seperated list of variables to define strata within which the variable is scaled", 
  nargs=Inf)
#### DESIGN MATRICES ####
p = add_argument(p, "--miDesignMatrix", help="components of design matrix for MI", nargs=Inf)
p = add_argument(p, "--seqDesignMatrix", help="components of design matrix for sequencing success", nargs=Inf)
#### MISSING VALUES ####
p = add_argument(p, "--miMissingDat", 
  help="list of data values to consider as missing. expected format is 
    var:missing_val;missing_min:missing_max:missing_prior1;missing_prior2;missing_std. Where 
    missing_val is the value to consider as missing, 
    missing_min is the minimum imputation value,
    missing_max is the maximum imputation value,
    missing_prior1 is the shape, 
    missing_prior2 is the scale, and 
    missing_std is the column to use to standardize missing values, use 'none' or omit for no standardization.", 
  nargs=Inf)
#### POST-HOC CALCULATIONS ####
p = add_argument(p, "--risks", nargs=Inf, help='post hoc calcualtion of risk')
p = add_argument(p, "--riskRatios", nargs=Inf,  help='post hoc calcualtion of risk ratios')
#### OUTPUT SETTINGS ####
p = add_argument(p, "--outAppend", help='additional string to append to out name', default='')
args = parse_args(p)

#args$dat = "output/211220_allreads_phsc_all_subgraphs_format_par.tsv"
#args$inputType = "full"
#args$stan = "stan/extended_model_hsp.stan"
#args$miDesignMatrix =c("round", "male_circumcision", "comm_type",
#  "sequencing_technology", "comm_type:sex", "comm_type:married", "comm_type:age_cat_coarse",
#  "comm_type:in_migrant", "comm_type:barworker")
#args$seqDesignMatrix =c("sequencing_technology", "sequencing_technology:log10_copies")
#args$scaleVars = 'log10_copies:sequencing_technology'
#args$risks = c("comm_type == 'inland'", "comm_type == 'fishing'")
#args$riskRatios =c("comm_type == 'fishing':comm_type == 'inland'")
#args$outAppend = "var_select"

#args$dat = "output/211220_allreads_phsc_all_subgraphs_format_par.tsv"
#args$inputType = "full"
#args$stan = "stan/full_model.stan"
#args$seqDesignMatrix =c("sequencing_technology", "sequencing_technology:log10_copies")
#args$scaleVars =c("log10_copies:sequencing_technology")
#args$risks =c("TRUE")
#args$outAppend = "test"

out = run_stan(args)
fit = out$fit
d_f = out$d_f %>% ungroup() %>% mutate(idx = seq(1,n()))
#fit = readRDS('fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_hsp_var_select.Rds')
#fit = readRDS('fit/extended_simulation_extended_model_hsp_.Rds')
# save input as Rds file

# summarise parameters
fit_draws = as_tibble(fit$draws(
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

fit_sum = bind_rows(
  fit_draws %>% 
    group_by(name) %>% 
    group_map(~get_hpd(.x$value) %>% 
      mutate(name=.y$name))) %>%
  left_join(fit_draws %>% group_by(name) %>%
    summarise(
      median = median(value),
      bulk_ess = ess_bulk(value),
      tail_ess = ess_tail(value),
      rhat = rhat(value)),
    by=c('name')) %>%
  select(name, median, lower, upper, bulk_ess, tail_ess, rhat)

id_cols = c('id', 'sample_id', 'study_id', 'round')
# add idx and add sample id to relevant parameters
fit_sum = fit_sum %>% mutate(
  idx = as.numeric(
    str_split(str_split(name, '\\[', simplify=TRUE)[,2], '\\]', simplify=TRUE)[,1])) %>%
  left_join(
    out$d_f %>% select(any_of(id_cols)) %>% mutate(idx = seq(1,n())), 
    by='idx')

if (!file.exists('fit/')){
  dir.create('fit') 
}

#### POST-HOC CALCULATIONS ####
#args$risks = c('logit_prob_MI+logit_prob_MI_coeffs[1]')
# add desired risks, if any
# get design matrix 
ind_prob_MI = as_tibble(fit$draws(
        'ind_log_prob_mi',
        inc_warmup = FALSE,
        format = "draws_df")) %>%
    mutate(cit = seq(1,n())) %>%
    pivot_longer(-cit) %>%
    filter(name != ".iteration" & 
      name != '.draw' & 
      name != '.chain') %>%
    mutate(idx = gsub("\\]", "", gsub("ind_log_prob_mi\\[", "", name))) %>%
    select(-name) %>%
    mutate(value = exp(value))

prob_MI = as_tibble(fit$draws(
      'prob_MI',
      inc_warmup = FALSE,
      format = "draws_df")) %>%
    mutate(cit = seq(1,n())) %>%
    pivot_longer(-cit) %>%
    filter(name != ".iteration" & 
      name != '.draw' & 
      name != '.chain') %>%
    mutate(idx = gsub("\\]", "", gsub("prob_MI\\[", "", name))) %>%
    select(-name)

if (!all(is.na(args$risks))){
  for (i in args$risks){
    i_ind_prob_MI = ind_prob_MI %>% 
      filter(idx %in% (d_f %>% filter(eval(parse(text=i))))$idx)
    i_prob_MI = prob_MI  %>% 
      filter(idx %in% (d_f %>% filter(eval(parse(text=i))))$idx)
    fit_sum = bind_rows(fit_sum, i_ind_prob_MI %>% 
        group_by(cit) %>%
        summarise(mean_value = mean(value)) %>%
        summarise(get_hpd(mean_value) %>%
            mutate(
              name = paste('sample_risk_', i, sep=''),
              median = median(mean_value),
              bulk_ess = ess_bulk(mean_value),
              tail_ess = ess_tail(mean_value),
              rhat = rhat(mean_value))))
    fit_sum = bind_rows(fit_sum, i_prob_MI %>% 
        group_by(cit) %>%
        summarise(mean_value = mean(value)) %>%
        summarise(get_hpd(mean_value) %>%
            mutate(
              name = paste('pop_risk_', i, sep=''),
              median = median(mean_value),
              bulk_ess = ess_bulk(mean_value),
              tail_ess = ess_tail(mean_value),
              rhat = rhat(mean_value))))
  }
}


if (!all(is.na(args$riskRatios))){
  for (i in args$riskRatios){
    i_split = str_split(i, ":", simplify=TRUE)
    num_ind_prob_MI = ind_prob_MI %>% 
      filter(idx %in% (d_f %>% filter(eval(parse(text=i_split[1]))))$idx) %>%
      group_by(cit) %>%
        summarise(mean_value = mean(value))
    num_prob_MI = prob_MI  %>% 
      filter(idx %in% (d_f %>% filter(eval(parse(text=i_split[1]))))$idx) %>%
       group_by(cit) %>%
        summarise(mean_value = mean(value))
    denom_ind_prob_MI = ind_prob_MI %>% 
      filter(idx %in% (d_f %>% filter(eval(parse(text=i_split[2]))))$idx) %>%
       group_by(cit) %>%
        summarise(mean_value = mean(value))
    denom_prob_MI = prob_MI  %>% 
      filter(idx %in% (d_f %>% filter(eval(parse(text=i_split[2]))))$idx) %>%
       group_by(cit) %>%
        summarise(mean_value = mean(value))
    fit_sum = bind_rows(
      fit_sum,
      num_ind_prob_MI %>% left_join(denom_ind_prob_MI, by=c('cit')) %>%
        mutate(value = mean_value.x / mean_value.y) %>%
        summarise(get_hpd(value) %>%
        mutate(
          name = paste('ind_rr_', i, sep=''),
          median = median(value),
          bulk_ess = ess_bulk(value),
          tail_ess = ess_tail(value),
          rhat = rhat(value))))
    fit_sum = bind_rows(
      fit_sum,
      num_prob_MI %>% left_join(denom_prob_MI, by=c('cit')) %>%
        mutate(value = mean_value.x / mean_value.y) %>%
        summarise(get_hpd(value) %>%
        mutate(
          name = paste('pop_rr_', i, sep=''),
          median = median(value),
          bulk_ess = ess_bulk(value),
          tail_ess = ess_tail(value),
          rhat = rhat(value))))
  }
}


# output file name
d_name = str_split(args$dat, '/', simplify=TRUE)
d_name = str_split(paste(d_name[2:length(d_name)], collapse='/'), '\\.', simplify=TRUE)
d_name = paste(d_name[1:length(d_name)-1], collapse='.')
s_name = str_split(str_split(args$stan, '/', simplify=TRUE)[2], '\\.', simplify=TRUE)
s_name = paste(s_name[1:length(s_name)-1], collapse='.')

op = paste(c('fit/', d_name, '_', s_name, '_', args$outAppend, '.Rds'), collapse='')
op_base = str_split(op, '/', simplify=TRUE)
op_base = paste(op_base[1:length(op_base)-1], collapse='/')

if (!dir.exists(op_base)){dir.create(op_base)}

# save output to rds
fit$save_object(file = op)
write_tsv(fit_sum, str_replace(op, ".Rds", "_summary.tsv"))
write_tsv(tibble(
  col=c(out$seq_design_matrix_cols, out$mi_design_matrix_cols), 
  type=c(
    rep('seq', length(out$seq_design_matrix_cols)),
    rep('mi', length(out$mi_design_matrix_cols)))),
  str_replace(op, ".Rds", "_design_cols.tsv"))
