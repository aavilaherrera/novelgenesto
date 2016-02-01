## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = '#>')

## ----init-chunk, message = FALSE, warning = FALSE------------------------
devtools::load_all('../package/novelgenesto/')

library(dplyr)
library(magrittr)
library(tidyr)
library(broom)
library(readr)

#library(ggplot2)
#library(GGally)
#library(gridExtra)

config_l = load_config('./taraoceans_config.json')
config_l

## ----input-chunk---------------------------------------------------------
# FUnkSFAM read counts
counts_df = load_counts(config_l$counts_fn)
#counts_df %<>% select(-starts_with('X'), starts_with('X') %>% sample(100))

presence_df = counts_to_presence(counts_df,
                                 threshold = config_l$counts_threshold
                                 )
# sample metadata
metadata_df = read_tsv(config_l$metadata_fn)

# add extra variables (fraction_cat, month, protocol_type)
metadata_df %<>% add_extra_vars()

data_df = inner_join(presence_df, metadata_df)

## ----clean-chunk---------------------------------------------------------
# drop extraneous columns
clean_df = data_df %>%
    select(run_accession,                      # sample identifier
           starts_with('X'),                   # FUnkSFAM identifiers
           one_of(unlist(config_l$predictor)), # environmental variables
           one_of(unlist(config_l$confounder)),        # confounders
           one_of(config_l$grouping)           # grouping variables
           ) %>%
    clean_99999(config_l$predictor$quantitative) %>%  # 99999 -> NA
    clean_99999(config_l$confounder$quantitative)

## ----qc-chunk------------------------------------------------------------
ff_stats_df = calc_FUnkSFAM_variation(clean_df, group = 'fraction_cat')
ph_stats_df = calc_phenotype_variation(clean_df, group = 'fraction_cat')
write_tsv(ff_stats_df, 'FUnkSFAM_variation.tsv')
write_tsv(ph_stats_df, 'phenotype_variation.tsv')

ff_stats_df
ph_stats_df

# subset most variable FUnkSFAMs for debuging
#ff_top_H = as.character(ff_stats_df$FUNKID) %>% unique() %>% head(10)
#clean_df = clean_df %>% select(-starts_with('X'), one_of(ff_top_H))
#clean_df

## ----model-chunk, warning = FALSE----------------------------------------
results_df = do_glm_tests(clean_df, 'fraction_cat', confounders = unlist(config_l$confounder))
results_df
write_tsv(results_df, 'results_raw.tsv')

## ------------------------------------------------------------------------
final_results_df = results_df %>% filter_results_all(ff_stats_df = ff_stats_df, config_l = config_l, group = 'fraction_cat') %>%
        adjust_pvalues() %>% format_final_results()
final_results_df
write_tsv(final_results_df, 'results_final.tsv')

