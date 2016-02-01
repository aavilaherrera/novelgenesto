# Functions to plot FUnkSFAM entropy

#' Gets 0\%, 25\%, 50\%, 75\%, 100\% quantiles of a column
#'
#' Note quantile and percentile used interchangeably here, oops.
#'
#' @param column A numeric column
#' @return A data.frame with percentiles and percentile values
#' @export
get_quantiles_df = function(column){
    quantile_df = data.frame(PercentileValue = quantile(column)) %>%
                        add_rownames(var = "Percentile")
    return(quantile_df)
}

#' Makes a histogram of FUnkSFAM entropies per group
#'
#' Uses ggplot2.
#'
#' @param ff_stats_df A data.frame of FUnkSFAM entropies
#'                     calculated by calc_FFvariation_per_group()
#' @param group Either 'HMP_BodySite' or 'HMP_Bodysubsite'
#' @return A ggplot2 object
#' @export
plot_FFentropy_by_group = function(ff_stats_df, group){
    quantile_df = ff_stats_df %>% group_by_(.dots = group) %>%
                        do({get_quantiles_df(.$Entropy_bits)})

    g = ggplot(ff_stats_df) + geom_histogram(aes(x = Entropy_bits)) +
        geom_vline(data = quantile_df,
                   aes(xintercept = PercentileValue,
                       color = Percentile,
                       linetype = Percentile
                       ),
                   size = 0.5,
                   show_guide = TRUE
                   ) +
        facet_wrap(as.formula(paste0('~', group))) +
        theme_bw() + theme(aspect.ratio = 1) +
        ggtitle("FUnkSFAM Presence-Absence Entropy")
    return(g)
}

#' @export
plot_scattermatrix = function(clean_df){
    clean_df %>% select(one_of(unlist(config_l$predictor),
                               config_l$covariate,
                               config_l$grouping
                               )
                        ) %>%
                 mutate_each(funs(as.factor), one_of(config_l$predictor$categorical,
                                                     config_l$grouping
                                                     )
                             ) %>%
                 plot(col = '#88888850')
}

#' @export
plot_ggpairs = function(df){
    if(ncol(df) > 5){
        warn("Too many columns for ggpairs()\n")
        colnames(df) = strtrim(colnames(df), 5)
    }
    gg = GGally::ggpairs(df)
    return(gg)
}

#' @export
plot_quantitative_predictors = function(df, group = NULL, quant_cols = config_l$predictor$quantitative){
    if(is.null(group)){
        melt_df = df %>% select(one_of(quant_cols)) %>%
            gather(Variable, Value, one_of(quant_cols))
        g = ggplot(melt_df, aes(x = Variable, y = Value))
    } else {
        melt_df = df %>% select(one_of(group, quant_cols)) %>%
            gather(Variable, Value, one_of(quant_cols), -one_of(group))
        g = ggplot(melt_df, aes_string(x = group, y = "Value"))
    }
    g = g + geom_violin(scale = 'width', trim = FALSE) +
            stat_summary(fun.y = median, geom = 'point', shape = 3) +
            facet_wrap(~Variable, scales = 'free')
    return(g)
}

#' @export
plot_categorical_predictors = function(df, group = NULL){
    cat_cols = config_l$predictor$categorical
    if(is.null(group)){
        melt_df = df %>% select(one_of(cat_cols)) %>%
            gather(Variable, Value, one_of(cat_cols))
        fgrp = '.'
    } else {
        melt_df = df %>% select(one_of(group, cat_cols)) %>%
            gather(Variable, Value, one_of(cat_cols), -one_of(group))
        fgrp = group
    }
    plt = function(mdf){
        g = ggplot(mdf, aes(x = Value)) + geom_bar() +
                facet_grid(as.formula(paste('Variable ~', fgrp)), scales = 'free', drop = TRUE) +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
                scale_x_discrete(labels = function(x){stringr::str_wrap(x, width = 40)})
        return(g)
    }
    glist = plyr::dlply(melt_df, plyr::.(Variable), plt)
    ml = do.call(gridExtra::marrangeGrob, args = c(glist, ncol = 1, nrow = 1))
    return(ml)
}

