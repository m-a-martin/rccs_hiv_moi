fit = readRDS('fit/211220_allreads_phsc_all_subgraphs_format_par_deep-phyloMI_var_select.Rds')
fit_draws = as_tibble(fit$draws(
	>     inc_warmup = FALSE,
	>     format = "draws_df"))
rowSums(fit_draws[,c("logit_prob_mi_coeffs[3]", "logit_prob_mi_coeffs[4]", "logit_prob_mi_coeffs[5]", "logit_prob_mi_coeffs[6]", "logit_prob_mi_coeffs[7]", "logit_prob_mi_coeffs[8]")])