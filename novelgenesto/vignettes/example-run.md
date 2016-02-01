The `novelgenesto` package accompanies supplemental methods for
**NOVELGENES\_TITLE\_HERE** by **AUTHORLIST**. It contains functions to
facilitate tests for association between environmental metadata from the
Tara Oceans expedition and under-annotated protein families investigated
in the manuscript (FUnkSFAMs).

See [novelgeneshmp](https://github.com/aavilaherrera/novelgeneshmp) for
a similar analysis in human microbiomes.

Finding associations between FUnkSFAM abundance and environmental metadata in Tara Oceans samples
-------------------------------------------------------------------------------------------------

Load some dependencies and helper packages...

    library(novelgenesto)

    library(dplyr)     # data wrangling
    library(magrittr)  #
    library(tidyr)     #
    library(broom)     # turn things into data.frames
    library(readr)     # read/write with sensible defaults

    # file handling for vignette
    vignette_datapath = function(path){
        extdata = system.file('extdata', package = 'novelgenesto')
        paste0(extdata, '/', basename(path))
    }

### Configuration file

The config file allows specification of predictor, outcome, and
confounding variables that are either categorical or quantitative. One
can also specify a grouping variable to stratify analyses (i.e. perform
the same analysis on non-overlapping groups of observations).

    # Specify configuration filename
    config_fn = vignette_datapath('taraoceans_configv.json')

    # Uses jsonlite::fromJSON to load the config into a list
    config_l = load_config(config_fn)

    config_l
    #> $counts_fn
    #> [1] "path/to/TO_counts.tsv"
    #> 
    #> $rpkg_fn
    #> [1] "path/to/TO_RPKG.tsv"
    #> 
    #> $metadata_fn
    #> [1] "path/to/tometadata_and_statistics.tsv.txt"
    #> 
    #> $span67_fn
    #> [1] "path/to/span_6_7_annotations.tsv"
    #> 
    #> $counts_threshold
    #> [1] 1
    #> 
    #> $predictor
    #> $predictor$quantitative
    #>  [1] "latitude_start"     "longitude_start"    "latitude_end"      
    #>  [4] "longitude_end"      "depth"              "temperature"       
    #>  [7] "salinity_sensor"    "oxygen_sensor"      "nitrate_sensor"    
    #> [10] "chlorophyll_sensor"
    #> 
    #> $predictor$categorical
    #> [1] "environment_biome"    "environment_feature"  "environment_material"
    #> [4] "protocol_type"       
    #> 
    #> 
    #> $grouping
    #> [1] "fraction_cat"
    #> 
    #> $confounder
    #> $confounder$quantitative
    #> [1] "count_reads"         "average_genome_size" "latitude_start"     
    #> 
    #> $confounder$categorical
    #> [1] "month"

Load the data
-------------

Load Tara Oceans metadata and FUnkSFAM presence-absence. The month a
sample was collected in (`month`) is derived from the sample collection
start date. `protocol_type` is parsed from `protocol_label`, and
`fraction_cat` manually groups size fractions with similar filter
thresholds (See function `add_extra_vars`).

    # vignette hack
    config_l$counts_fn = vignette_datapath(config_l$counts_fn)
    config_l$rpkg_fn = vignette_datapath(config_l$rpkg_fn)
    config_l$metadata_fn = vignette_datapath(config_l$metadata_fn)
    config_l$span67_fn = vignette_datapath(config_l$span67_fn)

    # FUnkSFAM read counts
    counts_df = load_counts(config_l$counts_fn)
    counts_df %<>% select(-starts_with('X'),
                          starts_with('X') %>% sample(100)  # subsample for vignette
                          )

    presence_df = counts_to_presence(counts_df,
                                     threshold = config_l$counts_threshold
                                     )
    # sample metadata
    metadata_df = read_tsv(config_l$metadata_fn)

    # add extra variables (fraction_cat, month, protocol_type)
    metadata_df %<>% add_extra_vars()

    data_df = inner_join(presence_df, metadata_df)
    #> Joining by: "run_accession"

Clean the data
--------------

    # drop extraneous columns
    clean_df = data_df %>%
        select(run_accession,                        # sample identifier
               starts_with('X'),                     # FUnkSFAM identifiers
               one_of(unlist(config_l$predictor)),   # environmental variables
               one_of(unlist(config_l$confounder)),  # confounders
               one_of(config_l$grouping)             # grouping variables
               ) %>%
        clean_99999(config_l$predictor$quantitative) %>%  # 99999 -> NA
        clean_99999(config_l$confounder$quantitative)

Quantify variation
------------------

### "Phenotype" (predictors, confounders) and FUnkSFAM presence-absence

The variation of FUnkSFAM presence-absence and other categorical
variables (e.g. `environment_biome`) is quantified with entropy.
Standard deviation is used for quantitative variables (e.g.
`nitrate_sensor`)

    ph_stats_df = calc_phenotype_variation(clean_df, group = 'fraction_cat')
    #> Warning: attributes are not identical across measure variables; they will
    #> be dropped
    #> Joining by: c("fraction_cat", "PHENONAME")
    #> Joining by: c("fraction_cat", "PHENONAME")
    #> Warning in outer_join_impl(x, y, by$x, by$y): joining factors with
    #> different levels, coercing to character vector
    ff_stats_df = calc_FUnkSFAM_variation(clean_df, group = 'fraction_cat')
    #> Warning: attributes are not identical across measure variables; they will
    #> be dropped
    #> Joining by: c("fraction_cat", "FUNKID")
    write_tsv(ph_stats_df, vignette_datapath('phenotype_variation.tsv'))
    write_tsv(ff_stats_df, vignette_datapath('FUnkSFAM_variation.tsv'))

    ph_stats_df
    #> Source: local data frame [68 x 5]
    #> 
    #>    fraction_cat           PHENONAME Entropy_bits ArbitraryStatistic
    #>           (chr)               (chr)        (dbl)              (lgl)
    #> 1       MultiGB         count_reads           NA                 NA
    #> 2        SmallG         count_reads           NA                 NA
    #> 3       MediumG         count_reads           NA                 NA
    #> 4          BigG         count_reads           NA                 NA
    #> 5        SmallG average_genome_size           NA                 NA
    #> 6          BigG average_genome_size           NA                 NA
    #> 7       MultiGB average_genome_size           NA                 NA
    #> 8       MediumG average_genome_size           NA                 NA
    #> 9          BigG               depth           NA                 NA
    #> 10      MediumG               depth           NA                 NA
    #> ..          ...                 ...          ...                ...
    #> Variables not shown: StandardDeviation (dbl)
    ff_stats_df
    #> Source: local data frame [400 x 4]
    #> 
    #>    fraction_cat    FUNKID Entropy_bits PercentPresent
    #>           (chr)    (fctr)        (dbl)          (dbl)
    #> 1       MediumG  X6_68255            1      0.6111111
    #> 2       MediumG  X5_58259            1      0.6111111
    #> 3        SmallG  X5_51930            1      0.5714286
    #> 4        SmallG  X5_15322            1      0.5714286
    #> 5          BigG X5_212159            1      0.5714286
    #> 6       MediumG  X5_31023            1      0.5555556
    #> 7       MediumG  X5_53861            1      0.5555556
    #> 8          BigG  X5_23582            1      0.5238095
    #> 9       MediumG  X6_64027            1      0.5000000
    #> 10         BigG  X5_59617            1      0.4761905
    #> ..          ...       ...          ...            ...

    # subset most variable FUnkSFAMs for debugging
    #ff_top_H = as.character(ff_stats_df$FUNKID) %>% unique() %>% head(10)
    #clean_df = clean_df %>% select(-starts_with('X'), one_of(ff_top_H))
    #clean_df

Fit logistic regression models and test for coefficients significantly larger than zero
---------------------------------------------------------------------------------------

Coefficients from a logistic regression were used to identify likely
associations between each environmental variable and presence of each
FUnkSFAM. The models account for latitude, month, read depth
(`count_reads`), and average genome size, and were fit for each group of
size fractions (`fraction_cat`) and for each pair of environmental
variable and FUnkSFAM.

### Test each predictor, accounting for confounders

This can take a long time, but is parallelizable (not implemented) as
all models are fit independently.

    results_df = do_glm_tests(clean_df, 'fraction_cat',
                              confounders = unlist(config_l$confounder)
                              )
    #> Joining by: "PHENONAME"
    #> Joining by: "PHENONAME"
    #> Joining by: "PHENONAME"
    #> Joining by: "PHENONAME"
    #> Joining by: "PHENONAME"
    #> Joining by: "PHENONAME"
    #> Joining by: "PHENONAME"
    #> Joining by: "PHENONAME"
    results_df
    #> Source: local data frame [47,400 x 21]
    #> 
    #>    fraction_cat                term      estimate    std.error
    #>           (chr)               (chr)         (dbl)        (dbl)
    #> 1          BigG         (Intercept) -2.556607e+01 7.544666e+05
    #> 2          BigG  chlorophyll_sensor -1.147681e-14 1.749481e+05
    #> 3          BigG         count_reads -2.392718e-23 1.832312e-03
    #> 4          BigG average_genome_size -1.167408e-21 6.281352e-02
    #> 5          BigG      latitude_start -2.887148e-17 1.976419e+04
    #> 6          BigG            monthJul  2.959142e-14 1.553917e+05
    #> 7          BigG            monthNov -3.327243e-15 4.313092e+05
    #> 8          BigG            monthOct  5.522660e-15 2.884739e+05
    #> 9          BigG            monthSep  1.471289e-16 3.792258e+05
    #> 10         BigG         (Intercept) -2.556607e+01 7.544666e+05
    #> ..          ...                 ...           ...          ...
    #> Variables not shown: statistic (dbl), p.value (dbl), null.deviance (dbl),
    #>   df.null (int), logLik (dbl), AIC (dbl), BIC (dbl), deviance (dbl),
    #>   df.residual (int), PHENONAME (chr), FUNKID (chr), N_samples (int),
    #>   FF_Entropy_bits (dbl), FF_PercentPresent (dbl), PH_Entropy_bits (dbl),
    #>   PH_ArbitraryStatistic (lgl), PH_StandardDeviation (dbl)
    write_tsv(results_df, vignette_datapath('results_raw.tsv'))

Filter the results and adjust p-values for false discovery
----------------------------------------------------------

### Filtering functions

<table>
<colgroup>
<col width="43%" />
<col width="56%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">name</th>
<th align="left">removes tests for...</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><code>filter_results_by_FFentropy_per_group</code></td>
<td align="left">mostly present or absent FUnkSFAMs</td>
</tr>
<tr class="even">
<td align="left"><code>filter_results_by_funksfam_annotation</code></td>
<td align="left">phylogenetically narrow FUnkSFAMs or with 10 or more annotations</td>
</tr>
<tr class="odd">
<td align="left"><code>filter_unneeded_tests</code></td>
<td align="left">confounder coefficients</td>
</tr>
<tr class="even">
<td align="left"><code>filter_results_all</code></td>
<td align="left">all of the above</td>
</tr>
</tbody>
</table>

    final_results_df = results_df %>%
            filter_results_all(ff_stats_df = ff_stats_df,
                               config_l = config_l,
                               group = 'fraction_cat'
                               ) %>%
            adjust_pvalues() %>%
            format_final_results()
    #> Joining by: "fraction_cat"
    #> Joining by: c("fraction_cat", "FUNKID")
    #> Warning in inner_join_impl(x, y, by$x, by$y): joining factor and character
    #> vector, coercing into character vector
    #> Joining by: "NoX_FFID"
    final_results_df
    #> Source: local data frame [60 x 18]
    #> 
    #>    fraction_cat   FUNKID           PHENONAME
    #>           (chr)    (chr)               (chr)
    #> 1       MultiGB X6_43108       oxygen_sensor
    #> 2       MultiGB X6_43108     longitude_start
    #> 3       MultiGB X6_43108  chlorophyll_sensor
    #> 4       MultiGB X6_43108      nitrate_sensor
    #> 5       MultiGB X6_43108       longitude_end
    #> 6       MultiGB X6_43108               depth
    #> 7       MultiGB X6_43108 environment_feature
    #> 8       MultiGB X6_43108         temperature
    #> 9       MultiGB X6_43108        latitude_end
    #> 10      MultiGB X6_43108 environment_feature
    #> ..          ...      ...                 ...
    #> Variables not shown: term (chr), p_adj_fdr (dbl), p.value (dbl), N_samples
    #>   (int), FF_Entropy_bits (dbl), FF_PercentPresent (dbl), PH_Entropy_bits
    #>   (dbl), PH_StandardDeviation (dbl), estimate (dbl), std.error (dbl),
    #>   statistic (dbl), logLik (dbl), df.residual (int), AIC (dbl), BIC (dbl)
    write_tsv(final_results_df, vignette_datapath('results_final.tsv'))
