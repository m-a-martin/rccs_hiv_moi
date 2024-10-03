suppressMessages(require(tidyverse))
suppressMessages(require(argparser))
suppressMessages(source('scripts/utils.R'))


p <- arg_parser("plot summary of empirical data")
p <- add_argument(p, "--dat", help="input data file", nargs=1)
p <- add_argument(p, "--fit", help="stan fit object", nargs=1)
p <- add_argument(p, "--colors", help="tsv file with color codes", default='config/colors.tsv', nargs=1)
p <- add_argument(p, "--out", help="output figure name", nargs=1)
args <- parse_args(p)
#args$fit = 'fit/211220_allreads_phsc_all_subgraphs_format_par_m_extended_model_sexpever_men.Rds'
#args$dat = 'output/211220_allreads_phsc_all_subgraphs_format_par_m.tsv'

cols = read_tsv(args$colors, show_col_types=FALSE)

args$colors_dict = setNames(cols$color, cols$var)

fit = readRDS(args$fit)

draws = as_tibble(fit$draws(
        inc_warmup = FALSE,
        format = "draws_df"))

# get just imputed sexpever values
imputed_sexpever = draws[grepl('X_mi_missing', colnames(draws)) & !grepl('raw', colnames(draws))] %>%
	mutate(cit = seq(1,n())) %>%
	pivot_longer(-cit) %>%
	mutate(idx = as.numeric(str_split(
		str_split(name, '\\[', simplify=TRUE)[,2],
		'\\]', simplify=TRUE)[,1]))

# read in summary so we can link with metadata
fit_sum = read_tsv(str_replace(args$fit, '.Rds', '_summary.tsv'), show_col_types=FALSE) 
idx_study_id = fit_sum %>% select(idx, study_id) %>% unique() %>% drop_na()

# read in input data 
dat = read_tsv(args$dat, show_col_types=FALSE)

dat = summarise_n_d(tabulate_n_d(dat %>% 
      filter(window_type == 'unique'))) %>%
	mutate(comm_type = ordered(comm_type, levels=c('inland', 'fishing'))) %>%
	left_join(idx_study_id, by='study_id') %>%
	filter(!is.na(idx)) %>%
	arrange(idx)

# read in design matrix
dm = read_tsv(str_replace(args$fit, '.Rds', '_design_cols.tsv'), show_col_types=FALSE) %>%
	filter(type == 'mi')

# replicate design matrix
# not! robust to other DMs
cols = unique(str_replace(str_replace(str_replace(dm$col, 'fishing', ''), 'inland', ''), 'bait_capture', ''))

mi_design_matrix = model.matrix(
      as.formula(
	paste('~', paste(cols, collapse='+'), collapse='')), data=dat)

where_missing = rbind(
	which(
        mi_design_matrix[,2:ncol(mi_design_matrix)] == 93,
        arr.ind=TRUE),
	which(
        mi_design_matrix[,2:ncol(mi_design_matrix)] == 100,
        arr.ind=TRUE))

missing_dat = dat[where_missing[,1],] %>%
	mutate(idx = seq(1,n()))

# merge in metadata
imputed_sexpever = imputed_sexpever %>% left_join(missing_dat, by='idx') %>%
	mutate(type = case_when(
		plhiv_sexpever_std == 93 ~ 'many',
		plhiv_sexpever_std == 100 ~ 'missing'))

# prior curves
groupsMany = dat %>% 
		select(plhiv_sexpever_std, comm_type, sex, age_cat_fine, 
			plhiv_sexpeverMany_shape, 
			plhiv_sexpeverMany_scale) %>%
		filter(plhiv_sexpever_std == 93) %>%
		unique()

xMany = seq(3,60,0.1)
fit_yMany = tibble(cbind(
	groupsMany[rep(seq(1,nrow(groupsMany)), each=length(xMany)),],
	tibble(x = rep(xMany, nrow(groupsMany))))) %>%
	# need to account for truncation
	mutate(y = 
		dlnorm(x, meanlog=plhiv_sexpeverMany_shape, sdlog=plhiv_sexpeverMany_scale) /
			(plnorm(60, meanlog=plhiv_sexpeverMany_shape, sdlog=plhiv_sexpeverMany_scale) - 
				plnorm(3, meanlog=plhiv_sexpeverMany_shape, sdlog=plhiv_sexpeverMany_scale)))

groupsMissing = dat %>% 
		select(plhiv_sexpever_std, comm_type, sex, age_cat_fine, 
			plhiv_sexpeverMissing_shape, 
			plhiv_sexpeverMissing_scale) %>%
		filter(plhiv_sexpever_std == 100) %>%
		unique()

xMissing = seq(1,60,0.1)
fit_yMissing = tibble(cbind(
	groupsMissing[rep(seq(1,nrow(groupsMissing)), each=length(xMissing)),],
	tibble(x = rep(xMissing, nrow(groupsMissing))))) %>%
	# need to account for truncation
	mutate(y = 
		dlnorm(x, meanlog=plhiv_sexpeverMissing_shape, sdlog=plhiv_sexpeverMissing_scale) /
			(plnorm(60, meanlog=plhiv_sexpeverMissing_shape, sdlog=plhiv_sexpeverMissing_scale) - 
				plnorm(1, meanlog=plhiv_sexpeverMissing_shape, sdlog=plhiv_sexpeverMissing_scale)))

fit_y = bind_rows(fit_yMany %>% mutate(type='many'), fit_yMissing %>% mutate(type = 'missing'))
p = ggplot(imputed_sexpever, aes(x=value, y=after_stat(density))) +
	geom_histogram(breaks=seq(0,61,1), fill='#eaeaea', color='#333333') +
	geom_line(data=fit_y, aes(x=x, y=y), color='#333333') +
	facet_grid(cols=vars(comm_type, type), rows=vars(age_cat_fine)) +
	xlab('lifetime sex partners') +
	ylab('density') +
	gtheme

ggsave(paste(c('figures/', args$out, '.pdf'), collapse=''), p, width=12, height=8)


