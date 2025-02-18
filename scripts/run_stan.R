suppressMessages(require(cmdstanr))
suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(require(HDInterval))
suppressMessages(require(posterior))
suppressMessages(source('scripts/utils.R'))



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
        print(d %>%
          group_by_at(group_vars_scale) %>% 
          summarise(m := mean(!!as.name(to_scale)), 
            s := sd(!!as.name(to_scale))))
        }else{
          d = d %>% 
            mutate(!!to_scale := scale(!!as.name(to_scale))[,1])
          print(d %>%
            summarise(m := mean(!!as.name(to_scale)), 
              s := sd(!!as.name(to_scale))))
        }
    }
  }
  return(d %>% ungroup())
}


format_data = function(args){
  # prepares input data for building stan data
  d = suppressMessages(read_tsv(args$dat) %>%
    filter(eval(parse(text=args$filter))))
  for (i in colnames(d)[grepl('round', colnames(d)) | grepl('risk_factor_MI', colnames(d))]){
    d[[i]] = as.character(d[[i]])
  }
  # tabulate and summarise data across windows
  if (args$inputType == 'full'){
    d = summarise_n_d(tabulate_n_d(d))
  }else if (args$inputType == "window"){
    d = summarise_n_d(d)
  }
  # scale variables that were requested to be scaled
  d = scale_vars(args$scaleVars, d)
  # remove rows with missing data
  vars = unique(Filter(function(x){x != "" & !is.na(x)}, 
    c(str_split(
          c(args$seqDesignMatrix,args$miDesignMatrix), "[:\\*]+", simplify=TRUE))))
  d_f = d[!apply(is.na(d[,vars]), 1, any),]
  n_rm = nrow(d)-nrow(d_f)
  if (n_rm == 0){n_rm="0"}else{n_rm=as.character(n_rm)}
  print(paste(n_rm, 
    ' rows dropped due to NA values', sep=''))
  return(d_f)
}


build_stan_model = function(args, stan_data){
  # pre-define shared functions
  build_data_param_blocks = function(label, stan_template, X, beta_cats){
    #### DATA BLOCK ####
    stan_template = str_replace(stan_template,
      gsub("LABEL", toupper(label), "// <!-- LABEL_MODEL_DATA -->"),
      gsub("LABEL", label, gsub("K", ncol(X), 
        "matrix[N_ind, K] X_LABEL; // LABEL input data ")))
    #### PARAMETERS BLOCK ####
    model_coeffs = ""
    beta_cats_counts = table(beta_cats)
    for (beta_cat in names(beta_cats_counts)){
      model_coeffs = paste0(model_coeffs,
          if_else(beta_cats_counts[beta_cat] == 1,
            paste0("  real LABEL_beta", beta_cat, ";\n"),
            paste0(gsub("K", beta_cats_counts[beta_cat], "  sum_to_zero_vector[K] LABEL_beta"),
              beta_cat, ";\n")))
    }
    # if shrinkage priors
    if (label == "mi" & args$miCoeffPriors == "shrinkage"){
      model_coeffs = paste0(model_coeffs,
          gsub("K", ncol(X),
            "\n  real<lower=0> r1_global;\n real<lower=0> r2_global;\n  vector<lower=0>[K] r1_local;\n  vector<lower=0>[K] r2_local;"))
    }
    stan_template = str_replace(stan_template,
      gsub("LABEL", toupper(label), "  // <!-- LABEL_MODEL_COEFFS -->"),
      gsub("LABEL", label, model_coeffs))
    #### TRANSFORMED-PARAMETERS BLOCK ####
    # add all coefficients into a single vector
    # if shrinkage priors
    if (label == "mi" & args$miCoeffPriors == "shrinkage"){
      # define vectors
      stan_template = str_replace(stan_template,
        gsub("LABEL", toupper(label), "  // <!-- LABEL_MODEL_COEFF_VEC -->"),
        gsub("LABEL", label, gsub("K", ncol(X),
          "  real<lower=0> tau;\n  vector<lower=0>[K] lambda;\n  vector[K] logit_prob_LABEL_coeffs;\n  tau = r1_global * sqrt(r2_global);\n  lambda = r1_local .* sqrt(r2_local);")))
      # build coeff vector
      n_per_cat = table(beta_cats)
      beta_labels = (tibble(cat = beta_cats, lag_cat = lag(cat)) %>%
        mutate(
          start = replace_na(cat != lag_cat, TRUE),
          id = cumsum(start)) %>%
        group_by(id) %>%
        mutate(x = paste0('[', row_number(), '],')))$x
      stan_template = str_replace(stan_template,
        gsub("LABEL", toupper(label), " // <!-- LABEL_MODEL_COEFF_VEC_ELEM -->"),
        gsub("LABEL", label, gsub(',\\]', '\\]', 
          paste0(
            '  // combine coefficients into a single vector\n  logit_prob_LABEL_coeffs = [',
            paste(paste0(paste0('LABEL_beta', beta_cats), 
            if_else(
              n_per_cat[beta_cats] > 1,
              beta_labels,
              ",")),
            collapse=''),
            "]' .* lambda*tau;"))))
    }else{
      # define coeff vector
      stan_template = str_replace(stan_template,
      gsub("LABEL", toupper(label), "  // <!-- LABEL_MODEL_COEFF_VEC -->"),
      gsub("LABEL", label, gsub("K", ncol(X), "  vector[K] logit_prob_LABEL_coeffs;")))
      # build coeff vector
      n_per_cat = table(beta_cats)
      beta_labels = (tibble(cat = beta_cats, lag_cat = lag(cat)) %>%
        mutate(
          start = replace_na(cat != lag_cat, TRUE),
          id = cumsum(start)) %>%
        group_by(id) %>%
        mutate(x = paste0('[', row_number(), '],')))$x
      stan_template = str_replace(stan_template,
        gsub("LABEL", toupper(label), " // <!-- LABEL_MODEL_COEFF_VEC_ELEM -->"),
        gsub("LABEL", label, gsub(',\\]', '\\]', 
          paste0(
            '  // combine coefficients into a single vector\n  logit_prob_LABEL_coeffs = [',
            paste(paste0(paste0('LABEL_beta', beta_cats), 
            if_else(
              n_per_cat[beta_cats] >1,
              beta_labels,
              ",")),
            collapse=''),
            "]';"))))
      }
      return(stan_template)
    }
  # read in template
  stan_template = readChar(args$stan, file.info(args$stan)$size)
  #### -------------------- ####
  #### IF SEQ SUCCESS MODEL ####
  #### -------------------- ####
  if (length(stan_data$X_seq) > 0){
    stan_template = 
      build_data_param_blocks("seq", stan_template, stan_data$X_seq, stan_data$seq_beta_cats)
    #### TRANSFORMED-PARAMETERS BLOCK ####
    # actual calculation of prob_seq
    stan_template = str_replace(stan_template,
      "  // <!-- SEQ_MODEL -->",
      paste0(
        "  prob_seq_1 = inv_logit( logit_prob_seq_baseline + logit_prob_seq_ind + X_seq*logit_prob_seq_coeffs );"))
    #### PRIORS BLOCK ####
    stan_template = str_replace(stan_template,
      "  // <!-- SEQ_MODEL_COEFF_PRIOR -->",
      "  target += normal_lpdf(logit_prob_seq_coeffs | 0, 2 );")
  }else{
    #### TRANSFORMED-PARAMETERS BLOCK ####
    # actual calculation of prob_seq
    stan_template = str_replace(stan_template,
      "  // <!-- SEQ_MODEL -->",
      paste0(
        "  prob_seq_1 = inv_logit( logit_prob_seq_baseline + logit_prob_seq_ind );"))
  }
  #### ------------ ####
  #### IF MI  MODEL ####
  #### ------------ ####
  if (length(stan_data$X_mi) > 0){
    stan_template = 
      build_data_param_blocks("mi", stan_template, stan_data$X_mi, stan_data$mi_beta_cats)
    #### TRANSFORMED-PARAMETERS BLOCK ####
    # actual calculation of prob_mi
    mi_model = paste0(
      "  matrix[N, K] X_mi_imputed = X_mi[,];\n",
      "  X_mi_missing = missing_min + (missing_max - missing_min) .* X_mi_missing_raw;\n",
      "  for (i in 1:N_missing){\n",
      "  X_mi_imputed[idx_missing[i,1], idx_missing[i,2]] = X_mi_missing[i] - X_mi_missing_std[i];\n",
      "  }\n",
      "  prob_mi = inv_logit( logit_prob_mi_baseline + X_mi_imputed*logit_prob_mi_coeffs );")
    stan_template = str_replace(stan_template,
      "  // <!-- MI_MODEL -->",
      gsub("K", ncol(stan_data$X_mi),
        gsub("N", nrow(stan_data$X_mi),
          gsub("N_missing", stan_data$N_missing,
            mi_model))))
    #### PRIORS BLOCK ####
    if (args$miCoeffPriors == "shrinkage"){
      nu = 2
      mi_priors = ''
      for (beta_cat in unique(stan_data$mi_beta_cats)){
        mi_priors = paste0(mi_priors, gsub("CAT", beta_cat, "  target += normal_lpdf(mi_betaCAT | 0, 1);\n"))
      }
      mi_priors = paste0(mi_priors,
        gsub("NU", nu,
          paste0("  // Wikipedia: the standard Cauchy distribution is the Student's t-distribution with one degree of freedom, and so it may be constructed by any method that constructs the Student's t-distribution.\n  // The location-scale t distribution results from compounding a Gaussian distribution (normal distribution) with mean and unknown variance, with an inverse gamma distribution placed over the variance with parameters a = nu/2 and b = nu*tau^2/2.\n",
            "  r1_local ~ normal(0.0, 1.0);\n  r2_local ~ inv_gamma(0.5*NU, 0.5*NU);\n  r1_global ~ normal(0.0, 1.0);\n  r2_global ~ inv_gamma(0.5, 0.5);")))
    }else{
      mi_priors = "  target += normal_lpdf(logit_prob_mi_coeffs | 0, 1 );"
    }
    stan_template = str_replace(stan_template,
      "  // <!-- MI_MODEL_COEFF_PRIOR -->",
      mi_priors)
  }else{
    #### TRANSFORMED-PARAMETERS BLOCK ####
    # actual calculation of prob_seq
    stan_template = str_replace(stan_template,
      "  // <!-- MI_MODEL -->",
      gsub("N", stan_data$N_ind,
        "  prob_mi = rep_vector( inv_logit( logit_prob_mi_baseline ), N );"))
  }
  return(stan_template)
}


format_missing_data = function(d, args, mi_design_matrix){
  # define output items
  missing = list()
  missing$N_missing = 0
  missing$idx_missing = matrix(0, 0,2)
  missing$X_mi_missing_std = vector('numeric')
  missing$missing_prior = vector('numeric')
  missing$missing_min = vector('numeric')
  missing$missing_max = vector('numeric')
  if (!all(is.na(args$miDesignMatrix))){
    if (!all(is.na(args$miMissingDat))){
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
          mi_design_matrix == missing_dat_val,
            MARGIN=2,
            grepl(missing_dat_col, colnames(mi_design_matrix)),
            '*') == 1, 
          arr.ind=TRUE)
        missing$idx_missing = rbind(missing$idx_missing, where_missing)
        # MINIMUM FOR IMPUTED VALUES
        missing$missing_min = c(
          missing$missing_min, rep(as.numeric(missing_data_min), nrow(where_missing)))
        # MAXIMUM FOR IMPUTED VALUES 
        missing$missing_max = c(
          missing$missing_max, rep(as.numeric(missing_data_max), nrow(where_missing)))
        # PRIOR FOR IMPUTED VALUES
        missing$missing_prior = rbind(missing$missing_prior, 
          as.matrix(d[where_missing[,1],c(missing_data_shape, missing_data_scale)]))
        # STANDARDIZATION FOR IMPUTED VALUES
        if (length(attr_split)){
          if (tolower(missing_data_std) != 'none' & !is.na(missing_data_std)){
            missing$X_mi_missing_std = c(
              missing$X_mi_missing_std, d[where_missing[,1],][[missing_data_std]])
          }else{
            missing$X_mi_missing_std = c(
              missing$X_mi_missing_std, rep(0, nrow(where_missing)))
          }
        }else{
          missing$X_mi_missing_std = c(
            missing$X_mi_missing_std, rep(0, nrow(where_missing)))}
        }
      missing$N_missing = nrow(missing$idx_missing)
    }
  }
  return(missing)
}


run_stan = function(args){
  #### DATA OBJECTS ####
  out = list()
  stan_data = list()
  #### READ IN AND PREPARE DATA ####
  d = format_data(args)
  out$d = d
  #### DATA ATTRIBUTES ####
  stan_data$N_ind = nrow(d)
  stan_data$N_obs_max = d$N_windows[1]
  stan_data$N_obs = d$N_obs
  stan_data$MI_obs = d$MI_obs
  #### BUILD SEQ DESIGN MATRIX ####
  seq_design_matrix = generate_dm(d, args$seqDesignMatrix)
  stan_data$X_seq = seq_design_matrix$dm
  if (length(seq_design_matrix$dm) > 0){stan_data$seq_beta_cats = seq_design_matrix$beta_cats}
  #### BUILD MI DESIGN MATRIX ####
  mi_design_matrix = generate_dm(d, args$miDesignMatrix)
  stan_data$nu = 2
  stan_data$X_mi = mi_design_matrix$dm
  if (length(mi_design_matrix$dm) > 0){stan_data$mi_beta_cats = mi_design_matrix$beta_cats}
  # handle missing data in MI predictors
  missing = format_missing_data(d, args, mi_design_matrix$dm)
  stan_data$N_missing = missing$N_missing
  stan_data$idx_missing = missing$idx_missing
  stan_data$X_mi_missing_std = missing$X_mi_missing_std
  stan_data$missing_prior = missing$missing_prior
  stan_data$missing_min = missing$missing_min
  stan_data$missing_max = missing$missing_max
  #### BUILD STAN MODEL TO MATCH DESIGN MATRICES ####
  stan_file = build_stan_model(args, stan_data)
  # save to file
  if (!dir.exists('stan/models')){dir.create('stan/models')}
  stan_file_path = gsub('.Rds', '.stan', gsub('fit/', 'stan/models/', args$op))
  writeLines(stan_file, stan_file_path)
  #### COMPILE STAN MODEL ####
  model_compiled = cmdstanr::cmdstan_model(stan_file_path)
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
    adapt_delta = if_else(length(mi_design_matrix$dm) == 0, 0.8, 0.95))
  out$stan_data = stan_data
  return(out)
}


#### --------------------- ####
#### MAIN CODE BEGINS HERE ####
#### --------------------- ####
p = arg_parser("run stan on input data")
#### INPUT DATA ####
p = add_argument(p, "--inputType", help="format of input file", default='summarized', type="character", nargs=1)
p = add_argument(p, "--dat", help="input data file", nargs=1)
p = add_argument(p, "--stan", help="stan model file", nargs=1)
p = add_argument(p, "--filter", default="TRUE", help="input data filter", nargs=1)
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
#### PRIORS ####
p = add_argument(p, "--miCoeffPriors", help="prior distributions for MI coeffs", default="normal")
#### POST-HOC CALCULATIONS ####
p = add_argument(p, "--strataPrevs", nargs=Inf, 
  help='post hoc calculation of risk based on design matrix')
p = add_argument(p, "--strataPrevRatios", nargs=Inf,  
  help='post hoc calculation of risk ratios based on design matrix')
p = add_argument(p, "--multivariateRisks", nargs=Inf, 
  help='post hoc calculation of risk assuming all individuals had a given set of covariates')
p = add_argument(p, "--multivariateRiskRatios", nargs=Inf,
  help='post hoc calculation of risk ratios assuming all individuals had a given set of covariates')
#### OUTPUT SETTINGS ####
p = add_argument(p, "--outAppend", help='additional string to append to out name', default='')
args = parse_args(p)
#args$dat = "simulations/extended_simulation.tsv"
#args$stan = "stan/deep-phyloMI.stan"
#args$miDesignMatrix =c("risk_factor_MI_1", "risk_factor_MI_2",    "risk_factor_MI_3", "risk_factor_MI_4", "risk_factor_MI_5")
#args$miCoeffPriors = "shrinkage"
#args$seqDesignMatrix = "scaled_log10_vl_obs"
#args$multivariateRiskRatios = "risk_factor_MI_1=='TRUE':risk_factor_MI_1=='FALSE'"


# output file name
d_name = str_split(args$dat, '/', simplify=TRUE)
d_name = str_split(paste(d_name[2:length(d_name)], collapse='/'), '\\.', simplify=TRUE)
d_name = paste(d_name[1:length(d_name)-1], collapse='.')
s_name = str_split(str_split(args$stan, '/', simplify=TRUE)[2], '\\.', simplify=TRUE)
s_name = paste(s_name[1:length(s_name)-1], collapse='.')

op = paste(c('fit/', d_name, '_', s_name, '_', args$outAppend, '.Rds'), collapse='')
args$op = op
out = run_stan(args)
fit = out$fit
d_f = out$d %>% ungroup() %>% mutate(idx = seq(1,n()))
#fit = readRDS('fit/211220_allreads_phsc_all_subgraphs_format_par_extended_model_hsp_var_select.Rds')
#fit = readRDS('fit/extended_simulation_extended_model_hsp_.Rds')
#fit = readRDS('fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_age_sex_comm.Rds')
# save input as Rds file
# summarise parameters
fit_draws = as_tibble(fit$draws(
    inc_warmup = FALSE,
    format = "draws_df"))

fit_draws = fit_draws[,
    colnames(fit_draws) %in%
      c('logit_prob_seq_baseline',
        'logit_prob_mi_baseline',
        'logit_prob_seq_ind_sd',
        'logit_prob_mi_fpr',
        'logit_prob_mi_fnr',
        'prob_mi_fpr',
        'prob_mi_fnr',
        'tau',
        'lambda') |
    grepl('logit_prob_seq_coeffs', colnames(fit_draws)) |
    grepl('logit_prob_mi_coeffs', colnames(fit_draws)) |
    grepl('lambda', colnames(fit_draws)) |
    (grepl('prob_mi\\[', colnames(fit_draws)) & 
      !grepl('log', colnames(fit_draws)))]  %>%
  mutate(cit = seq(1,n()))


fit_sum = fit_draws %>% select(-cit) %>% reframe(across(all_of(colnames(.)), 
    ~c(hdi(., credMass=0.95), 
      median(.),
      ess_bulk(.),
      ess_tail(.),
      rhat(.)))) %>%
  mutate(i = c('lower', 'upper', 'median', 'bulk_ess', 'tail_ess', 'rhat')) %>%
  pivot_longer(-i) %>%
  pivot_wider(names_from=i, values_from=value)


summarize_group = function(x, name){
  return(get_hpd(x$value) %>%
              mutate(
                name = name,
                median = median(x$value),
                bulk_ess = ess_bulk(x$value),
                tail_ess = ess_tail(x$value),
                rhat = rhat(x$value)))
}

id_cols = c('id', 'sample_id', 'study_id', 'round')
# add idx and add sample id to relevant parameters
fit_sum = fit_sum %>% mutate(
    idx = if_else(
    grepl('prob_mi\\[', name),
    str_split(str_split(name, '\\[', simplify=TRUE)[,2], '\\]', simplify=TRUE)[,1],
    NA)) %>%
  left_join(d_f %>% mutate(idx = as.character(idx)) %>% select(any_of(c(id_cols, 'idx'))),
    by='idx')


#### POST-HOC CALCULATIONS ####
# for each unique strata defined by the design matrix calculate the probability of multiple infection
# first get unique data rows based on design matrix and get Ns 
uniq_dat = d_f %>%
  group_by_at(unique(Filter(function(x){x != ""}, c(str_split(args$miDesignMatrix, "\\*", simplify=TRUE))))) %>%
  mutate(n=n()) %>%
  slice(1) %>%
  ungroup() %>%
  select(any_of(c('idx', unique(Filter(function(x){x != ""}, c(str_split(args$miDesignMatrix, "\\*", simplify=TRUE)))), 'n')))
# then get prob_MI for each iteration associated with those unique data rows
uniq_prob_mi = bind_cols(
  fit_draws[,'cit'],
  fit_draws[
    grepl('prob_mi\\[', colnames(fit_draws)) & !grepl('logit', colnames(fit_draws))][
    ,uniq_dat$idx])

# calculate average prob_MI in the population for each iteration
prob_mi = uniq_prob_mi %>%
  pivot_longer(-cit) %>%
  mutate(idx = as.numeric(gsub("\\]", "", gsub("prob_mi\\[", "", name)))) %>%
  left_join(uniq_dat, by='idx')

# add average MI to summary
fit_sum = bind_rows(
  fit_sum,
  prob_mi %>% 
  group_by(cit) %>% 
  summarise(value = sum(value*n)/sum(n)) %>%
  summarise(summarize_group(., 'prob_mi')))

# strata risks
if (!all(is.na(args$strataPrevs))){
  for (i in args$strataPrevs){
    fit_sum = bind_rows(
      fit_sum,
      prob_mi %>%
        filter(eval(parse(text=i))) %>%  
        group_by(cit) %>%
        summarise(value = sum(value*n)/sum(n)) %>%
        summarise(summarize_group(., paste(c('prob_mi_', gsub('"', "", i)), collapse=''))))
  }
}

# strata risk ratios
if (!all(is.na(args$strataPrevRatios))){
  for (i in args$strataPrevRatios){
    i_split = str_split(i, ":", simplify=TRUE)
    fit_sum = bind_rows(
      fit_sum,
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
        summarise(summarize_group(., paste(c('pr_mi', "_", gsub('"', "", i)), collapse=''))))
  }
}

# multivariate risks
if (!all(is.na(args$multivariateRisks))){
  for (i in args$multivariateRisks){
    # get intercept and coefficients
    coeffs = get_coeffs(fit_draws)
    k_split = str_split(i, "&", simplify=TRUE)
    # pull original unique data frame
    k_df = uniq_dat
    for (j in k_split){
      j_split = str_split(j, "==", simplify=TRUE)
      k_df[gsub("[[:blank:]]", "", j_split[1])] = eval(parse(text=j_split[2]))
    }
    # bind dataframes as a hack to get correct design matrices when only one level
    dm = generate_dm(bind_rows(k_df,uniq_dat), args$miDesignMatrix)$dm[1:nrow(k_df),]
    # calculate risk
    risk = apply(as.matrix(coeffs)[,2:ncol(coeffs)], 1, 
      function(x) inv_logit(colSums(t(cbind(1, dm))*x)))
    mean_risk = as_tibble(colSums(risk*uniq_dat$n)/sum(uniq_dat$n)) %>%
      mutate(cit = seq(1,n()))
    fit_sum = bind_rows(
      fit_sum,
      mean_risk %>%
        summarise(summarize_group(., paste('multivar_risk_mi_', i, sep=''))))
  }
}


if (!all(is.na(args$multivariateRiskRatios))){
  for (i in args$multivariateRiskRatios){
    i_split = str_split(i, ":", simplify=TRUE)
    # get intercept and coefficients
    coeffs = get_coeffs(fit_draws)
    # create design matrices for desired RR, keeping all other covariates the same
    # risk ratio so only should be two items in i
    dms = list()
    for (k in 1:2){
      k_split = str_split(i_split[k], "&", simplify=TRUE)
      # pull original data frame
      k_df = uniq_dat
      for (j in k_split){
        j_split = str_split(j, "==", simplify=TRUE)
        k_df[gsub("[[:blank:]]", "", j_split[1])] = eval(parse(text=j_split[2]))
      }
      # bind dataframes as a hack to get correct design matrices when only one level
      dms[[k]] = generate_dm(bind_rows(k_df,uniq_dat), args$miDesignMatrix)$dm[1:nrow(k_df),]
    }
    # calculate risk for each dm
    risks = tibble()
    for (k in 1:length(dms)){
      k_risk = apply(as.matrix(coeffs)[,2:ncol(coeffs)], 1, 
        function(x) inv_logit(colSums(t(cbind(1, dms[[k]]))*x)))
      risks = bind_rows(risks,
        as_tibble(colSums(k_risk*uniq_dat$n)/sum(uniq_dat$n)) %>%
          mutate(
            k=k,
            cit = seq(1,n())))
    }
    # for each iteration calculate RR
    rr = risks %>% pivot_wider(names_from=k, values_from=value) %>%
      mutate(value = `1`/`2`)
    fit_sum = bind_rows(
      fit_sum,
      rr %>% summarise(summarize_group(., paste('multivar_rr_mi_', i, sep=''))))
  }
}


if (!file.exists('fit/')){
  dir.create('fit') 
}


op_base = str_split(op, '/', simplify=TRUE)
op_base = paste(op_base[1:length(op_base)-1], collapse='/')

if (!dir.exists(op_base)){dir.create(op_base)}

# save output to rds
fit$save_object(file = op)
write_tsv(fit_sum, str_replace(op, ".Rds", "_summary.tsv"))
design_cols = tibble()


if (!any(is.na(args$seqDesignMatrix))){
  seq_cols = colnames(out$stan_dat$X_seq)
  design_cols = bind_rows(design_cols,
    tibble(
      col=seq_cols,
      type=rep('seq', length(seq_cols))))
}
if (!any(is.na(args$miDesignMatrix))){
  mi_cols = colnames(out$stan_dat$X_mi)
  design_cols = bind_rows(design_cols,
    tibble(
      col=mi_cols,
      type=rep('mi', length(mi_cols))))
}




write_tsv(design_cols,
  str_replace(op, ".Rds", "_design_cols.tsv"))


