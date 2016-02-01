# Function to load and save various files and file formats

#' Load abundance data (as counts)
#'
#' Loads and transposes counts csv
#' where each row is a FUnkSFAM. The first column is a FUNKSFAM ID,
#' other columns are counts for each run_accession
#'
#' @param fn A filename
#' @return A data.frame with each row a sample and FUnkSFAMs as columns
#' @export
load_counts = function(fn){
    #cnts_df = load_csv(fn)
    cnts_df = readr::read_tsv(fn)  # Tara Oceans counts are tab separated
    colnames(cnts_df)[1] = 'FUNKID'
    cnts_df$X.1 = NULL  # There was an extra column to delete
    cnts_df %<>% gather(key = 'run_accession', value = 'Counts', -FUNKID) %>%
              select(FUNKID, run_accession, Counts) %>%
              mutate(FUNKID = paste('X', FUNKID, sep = '')) %>%  # "X" preprended to colnames
              spread(FUNKID, Counts) %>%
              mutate(run_accession = as.character(run_accession)) %>%    # Enforce types
              mutate_each(funs(as.integer), starts_with('X'))  # Enforce types
    return(cnts_df)
}

#' Formats final results
#'
#' Selects interesting results to print, sorts by adjusted and
#' unadjusted p-values
#'
#' @param res_df A data.frame of results
#' @return A formatted res_df
#' @export
format_final_results = function(res_df, GRP = 'fraction_cat'){
    # choose columns and order to report in table
    res_df %<>% select(matches(GRP),
                       FUNKID, PHENONAME, term, p_adj_fdr, p.value,
                       N_samples,
                       FF_Entropy_bits, FF_PercentPresent,
                       PH_Entropy_bits,
                       PH_StandardDeviation,
                       estimate, std.error, statistic,
                       logLik, df.residual, AIC, BIC
                       ) %>%
                arrange(p_adj_fdr, p.value)
    return(res_df)
}
