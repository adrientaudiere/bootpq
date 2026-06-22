# bootpq 0.0.0 (development version)

* Initial development version of the package, providing bootstrapping by rarefaction permutations for 'phyloseq' objects.
* Add a "Get started with bootpq" vignette and a pkgdown website skeleton.
* `adonis_rarperm_pq()` computes PERMANOVA (via `MiscMetabar::adonis_pq()`) across many rarefaction permutations and summarizes the results as mean and quantiles, migrated from 'MiscMetabar'.
* `hill_test_rarperm_pq()` tests the effect of a factor on Hill numbers across many rarefaction permutations using `ggstatsplot::ggbetweenstats()`, migrated from 'MiscMetabar'.
* `var_par_rarperm_pq()` partitions the variation of a 'phyloseq' object across many rarefaction permutations, averaging adjusted R squared and reporting quantiles, migrated from 'MiscMetabar'.
