# Functions for cleaning metadata

#' Sets 99999 -> NA
#'
#' @export
clean_99999 = function(meta_df, dirty_cols){
    clean_meta_df = meta_df %>% mutate_each(funs({ifelse(. == 99999, NA, .)}),
                                            one_of(dirty_cols)
                                            )
    return(clean_meta_df)
}

