# Functions for loading and cleaning subject metadata

#' Add extra variables
#'
#' Combines fraction size ranges into categories `fraction_cat`,
#' assigns a protocol type by parsing the first token
#' of `protocol_label`, and extracts the month from `event_date_time_start`.
#'
#' @param meta_df A data.frame of Tara Oceans metadata
#' @return meta_df with extra categories
#' @export
#'
add_extra_vars = function(meta_df){
    meta_df %<>% mutate(fraction_cat = ifelse(size_fraction_lower_threshold == 0.1 &
                                             size_fraction_upper_threshold == 0.22,
                                                'SmallG',
                                        ifelse(size_fraction_lower_threshold == 0.22 &
                                               size_fraction_upper_threshold == 0.45,
                                                'MediumG',
                                        ifelse(size_fraction_lower_threshold == 0.45 &
                                               size_fraction_upper_threshold == 0.8,
                                                'BigG',
                                        ifelse(size_fraction_lower_threshold == 0.22 &
                                               size_fraction_upper_threshold %in% c(1.6, 3),
                                                'MultiGB',
                                                'Other'
                                        ))))
                        )
   meta_df %<>% separate(protocol_label, 'protocol_type', sep = '_', extra = 'drop')
   meta_df %<>% mutate(month = as.character(lubridate::month(as.POSIXct(event_date_time_start), label = TRUE)))
   return(meta_df)
}
