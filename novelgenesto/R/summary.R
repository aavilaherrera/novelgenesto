# Functions for computing FUnkSFAM and phenotype stats

#' Gets percent TRUE in a presence-absence column
#'
#' @param ff_pres A presence-absence logical vector
#' @return A percentage
#' @export
get_percent_present = function(ff_pres){
    return(sum(ff_pres) / length(ff_pres))
}

#' Gets entropy (James-Stein shrinkage estimator) for a categorical column
#'
#' Requires \code{entropy} package
#'
#' @param cat_column A categorical vector
#' @return Entropy in bits
#' @export
get_entropy = function(cat_column){
    return(entropy::entropy(table(cat_column),
                            verbose = FALSE,
                            method = 'shrink',
                            unit = 'log2'
                            )
           )
}

#' Gets ArbitraryStatistic for a categorical column
#'
#' I don't know what to call this.
#'
#' Checks that there is more than 1 value with more than 4 counts
#'
#' @param cat_column A categorical vector
#' @return ArbitraryStatistic
#' @export
get_arbitrary_statistic = function(cat_column){
    return(sum(table(cat_column) > 4) > 1)
}

#' Get number of samples per group (body site or subsite)
#'
#' @param df A data.frame in `tidy` format
#' @param group Either 'HMP_BodySite' or 'HMP_BodySubsite'
#' @return A data.frame
#' @export
calc_nsamples_per_group = function(df, group){
    nsamp_df = df %>% group_by_(.dots = group) %>% tally()
    return(nsamp_df)
}

#' Get FUnkSFAM presence-absence variation per group
#'
#' @inheritParams calc_nsamples_per_group
#' @return A data.frame wih Entropy_bits and PercentPresent per FUNKID per group
#' @export
calc_FUnkSFAM_variation = function(df, group = NULL){
    df %<>% group_by_(.dots = group)
    entropy_df = df %>% summarise_each(funs(get_entropy), starts_with('X')) %>%
                        gather(FUNKID, Entropy_bits, starts_with('X'))
    percent_pres_df = df %>% summarise_each(funs(get_percent_present), starts_with('X')) %>%
                        gather(FUNKID, PercentPresent, starts_with('X'))
    variation_df = full_join(entropy_df, percent_pres_df) %>%
                        arrange(desc(Entropy_bits), desc(PercentPresent))
    return(variation_df)
}

#' Get phenotype variation per group
#'
#' @inheritParams calc_nsamples_per_group
#' @return A data.frame with standard deviation, Entropy_bits, and ArbitraryStatistic per PHENONAME
#'         per group. NAs if variation statistic is not appropriate.
#' @export
calc_phenotype_variation = function(df, group = NULL){
    df %<>% group_by_(.dots = group)
    quantitative = unlist(c(config_l$predictor$quantitative, config_l$confounder$quantitative))
    categorical = unlist(c(config_l$predictor$categorical, config_l$confounder$categorical))

    # compute quantitative variation statistics
    quant_cols = base::intersect(colnames(df), quantitative)
    quant_var_df = df %>% summarise_each(funs({sd(., na.rm = TRUE)}), one_of(quant_cols)) %>%
                          gather(PHENONAME, StandardDeviation, one_of(quant_cols))

    cat_cols = base::intersect(colnames(df), categorical)
    entropy_df = df %>% summarise_each(funs(get_entropy), one_of(cat_cols))
    arbitrary_stat_df = df %>% summarise_each(funs(get_arbitrary_statistic), one_of(cat_cols))

    cat_var_df = full_join(gather(entropy_df, PHENONAME, Entropy_bits,
                                  one_of(cat_cols)
                                  ),
                           gather(arbitrary_stat_df, PHENONAME, ArbitraryStatistic,
                                  one_of(cat_cols)
                                  )
                           )
    variation_df = full_join(cat_var_df,
                             quant_var_df
                             ) %>%
                        arrange(desc(StandardDeviation),
                                desc(Entropy_bits), desc(ArbitraryStatistic)
                                )
    return(variation_df)
}
