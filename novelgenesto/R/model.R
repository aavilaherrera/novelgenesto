# Functions for logistic regression modeling and testing coefficients

#' Make all pairs of FUNKID and PHENONAME from variable column names
#'
#' The column names of a data.frame in `tidy` format contain
#' many FUnkSFAM IDs (FUNKIDs) and many phenotypes variables (PHENONAME)
#' along with other metadata. This function extracts pairs of FUNKID
#' and PHENONAME.
#'
#' Assumes only FUNKIDs start with 'X'. Assumes \code{config_l} exists
#' @param column_names A character vector
#' @return A data.frame of (FUNKID, PHENONAME) pairs.
#' @export
pick_tests = function(column_names){
    ff_names = grep('^X', column_names, value = TRUE)
    ph_names = base::intersect(column_names, unlist(config_l$predictor))
    ff_ph = setNames(expand.grid(ff_names, ph_names),
                     c('FUNKID', 'PHENONAME')
                     ) %>%
            mutate_each(funs(as.character))

    #warn("DBG: SUBSAMPLING FOR TESTING PURPOSES\n")
    #if(nrow(ff_ph) > 20){
    #    ff_ph %<>% sample_n(20)
    #}

    return(ff_ph)
}

#' Fits glm model, gets coefficient stats, and formats results
#'
#' NOTE: if glm() fails, this test is silently skipped
#'
#' @param df A data.frame in `tidy` format (rows are observations, columns are
#'           variables)
#' @param ff_name FUNKID name to test
#' @param ph_name PHENONAME to test
#' @return A data.frame with statistics on coefficients, input variation, and
#'         sample size
#' @export
fit_logistic = function(df, ff_name, ph_name, confounders = NULL){
    rhs = paste(c(ph_name, confounders), collapse = '+')
    the_formula = as.formula(paste(ff_name, '~', rhs))
    the_model = try(
                    glm(formula = the_formula, data = df,
                        family = 'binomial'
                        ),
                    silent = TRUE  # suppresses error messages
                    )
    if('try-error' %ni% class(the_model)){
        coeff_df = tidy(the_model)
        gof_df = glance(the_model)
        stats_df = summarize_model(the_model, ff_name, ph_name)
        res_df = Reduce(merge, list(coeff_df, gof_df, stats_df)) %>% tbl_df()
        return(res_df)
    } else {
        warn(paste("WW: There was an error testing:",
                   ff_name, ph_name, "\n",
                   collapse = " "
                   )
             )
        return(data.frame())
    }
}

#' Runs \link{\code{fit_logistic}} for each (FUNKID, PHENONAME) pair.
#'
#' @param df A data.frame in `tidy` format (rows are observations, columns are
#'           variables)
#' @return A data.frame with statistics on coefficients, input variation, and
#'         sample size
#' @export
glm_loop = function(df, confounders = NULL){
    df %<>% prefilter_by_arbitrary_statistic(calc_phenotype_variation(df))
    ff_ph = pick_tests(colnames(df))
    warn(paste0("\nAbout to fit ", nrow(ff_ph), " models!\n"))
    res_df = ff_ph %>% rowwise() %>%
                    do({fit_logistic(df, .$FUNKID, .$PHENONAME, confounders)})
    return(res_df)
}

#' Tests for significant associations between FUNKIDs and PHENONAMEs per grouping
#'
#' Runs \link{\code{glm_loop}} for each \code{HMP_BodySite} or \code{HMP_BodySubsite}
#'
#' @param df A data.frame in `tidy` format (rows are observations, columns are
#'           variables)
#' @return A data.frame with statistics on coefficients, input variation, and
#'         sample size
#' @export
do_glm_tests = function(df, group_name, confounders = NULL){
    res_df = df %>% group_by_(.dots = group_name) %>%
                do({glm_loop(., confounders)}) %>%
                ungroup()
    return(res_df)
}

#' Get statistics from the cleaned model data.frame in glm()
#'
#' Gets Entropy (bits) and ArbitraryStatistic for \code{ph_name}
#' Gets Entropy (bits) and PercentPresent for \code{ff_name}
#' Get N_samples (number of rows in glm()$model data.frame)
#'
#' DBG NOTE: \code{ff_name} and \code{ph_name} could probably be extracted from
#'           \code{names(the_model)}
#'
#' @param the_model The data.frame from glm()$model
#' @param ff_name FUNKID tested
#' @param ph_name PHENONAME tested
#' @return A data.frame with the model statistics
#' @export
summarize_model = function(the_model, ff_name, ph_name){
    ph_values = the_model$model[, ph_name]
    if(is.numeric(ph_values)){
        PH_Entropy_bits = NA
        PH_ArbitraryStatistic = NA
        PH_StandardDeviation = sd(ph_values, na.rm = TRUE)
    } else {
        PH_Entropy_bits = ph_values %>% get_entropy()
        PH_ArbitraryStatistic = ph_values %>% get_arbitrary_statistic()
        PH_StandardDeviation = NA
    }
    stats_df = data.frame(PHENONAME = ph_name, FUNKID = ff_name,
                          N_samples = the_model$model %>% nrow(),
                          FF_Entropy_bits = the_model$model[, ff_name] %>% get_entropy(),
                          FF_PercentPresent = the_model$model[, ff_name] %>% get_percent_present(),
                          PH_Entropy_bits = PH_Entropy_bits,
                          PH_ArbitraryStatistic = PH_ArbitraryStatistic,
                          PH_StandardDeviation = PH_StandardDeviation
                          )
    return(stats_df)
}

