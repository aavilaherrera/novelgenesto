## About

This package accompanies the manuscript NOVELGENES_TITLE_HERE. It includes
functions used to clean, filter, and analyze data for associations between
FUnkSFAM presence and subject metadata in the manuscript. See the
[vignette](./novelgenesto/vignettes/example-run.md) for an example.

## Install

Download the tarball and install with

```bash

R CMD INSTALL novelgenesto_0.0.0.9000.tar.gz

```

Or with [devtools](https://github.com/hadley/devtools)

```r

library(devtools)
install_github('aavilaherrera/novelgenesto@dev', subdir = 'novelgenesto')

# To build vignettes locally (takes a few extra seconds):
install_github('aavilaherrera/novelgenesto@dev',
               subdir = 'novelgenesto',
               build_vignettes = TRUE
               )

```

## Results Summary

Results for the whole data set are summarized in
[results_summary](./results_summary/) as the raw results are ~500MB and ~2
million lines long.

## License

GPLv3. See [LICENSE](./novelgenesto/LICENSE).
