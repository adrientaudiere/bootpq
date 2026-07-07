# Get started with bootpq

## Why bootpq?

A single rarefaction to even sequencing depth throws away part of the
data and depends on one random draw. The conclusion of a downstream test
(PERMANOVA, a Hill-number group comparison, a variation partitioning)
can therefore change from one rarefaction to the next.

`bootpq` removes this fragility by **repeating the analysis across many
rarefaction permutations** — each with a different random seed — and
summarising the results (mean and quantiles). You get an estimate
together with its variability across rarefactions, instead of a single,
possibly atypical, value.

`bootpq` is part of the [pqverse](https://github.com/adrientaudiere)
ecosystem and builds on
[`MiscMetabar`](https://adrientaudiere.github.io/MiscMetabar/).

## Installation

``` r
# install.packages("remotes")
remotes::install_github("adrientaudiere/bootpq")
```

## A worked example

We use `data_fungi_mini`, a small fungal metabarcoding dataset shipped
with `MiscMetabar`. We first drop samples with missing covariates.

``` r
library(bootpq)
library(MiscMetabar)
data(data_fungi_mini)

data_fungi_woNA <- subset_samples(
  data_fungi_mini,
  !is.na(Time) & !is.na(Height)
)
data_fungi_woNA
#> phyloseq-class experiment-level object
#> otu_table()   OTU Table:         [ 45 taxa and 75 samples ]
#> sample_data() Sample Data:       [ 75 samples by 7 sample variables ]
#> tax_table()   Taxonomy Table:    [ 45 taxa by 12 taxonomic ranks ]
#> refseq()      DNAStringSet:      [ 45 reference sequences ]
```

### PERMANOVA across rarefaction permutations

[`adonis_rarperm_pq()`](https://adrientaudiere.github.io/bootpq/reference/adonis_rarperm_pq.md)
runs
[`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html)
on each rarefied version of the object and summarises the results. We
keep `nperm` small here so the vignette builds quickly; in practice use
the default (`nperm = 99`) or more.

``` r
res <- adonis_rarperm_pq(
  data_fungi_woNA,
  "Time*Height",
  na_remove = TRUE,
  nperm = 9,
  progress_bar = FALSE
)
#> Taxa are now in columns.
#> Removing NA from Time
#> Removing NA from Height
#> Taxa are now in columns.
#> Removing NA from Time
#> Removing NA from Height
#> Taxa are now in columns.
#> Removing NA from Time
#> Removing NA from Height
#> Taxa are now in columns.
#> Removing NA from Time
#> Removing NA from Height
#> Taxa are now in columns.
#> Removing NA from Time
#> Removing NA from Height
#> Taxa are now in columns.
#> Removing NA from Time
#> Removing NA from Height
#> Taxa are now in columns.
#> Removing NA from Time
#> Removing NA from Height
#> Taxa are now in columns.
#> Removing NA from Time
#> Removing NA from Height
#> Taxa are now in columns.
#> Removing NA from Time
#> Removing NA from Height

# Mean of the PERMANOVA table across the 9 rarefaction permutations
res$mean
#>          Df  SumOfSqs        R2        F    Pr(>F)
#> Model     5  2.429663 0.0687147 1.019227 0.4693333
#> Residual 69 32.933300 0.9312853       NA        NA
#> Total    74 35.362963 1.0000000       NA        NA
```

`res$quantile_min` and `res$quantile_max` give the corresponding lower
and upper quantiles, so you can see how stable each term is across
rarefactions.

### Hill-number group tests

[`hill_test_rarperm_pq()`](https://adrientaudiere.github.io/bootpq/reference/hill_test_rarperm_pq.md)
compares Hill numbers (q = 0, 1, 2) between the levels of a factor
across rarefaction permutations. It relies on `ggstatsplot`, so the
chunk below is shown but not run during the vignette build.

``` r
hill_test_rarperm_pq(data_fungi_woNA, fact = "Height", nperm = 99)
```

### Variation partitioning

[`var_par_rarperm_pq()`](https://adrientaudiere.github.io/bootpq/reference/var_par_rarperm_pq.md)
repeats a variation partitioning across rarefaction permutations,
summarising how much variation each set of explanatory variables
explains.

``` r
var_par_rarperm_pq(
  data_fungi_woNA,
  list_component = list(
    "Time" = c("Time"),
    "Size" = c("Height", "Diameter")
  ),
  nperm = 99
)
```

## Where to go next

- The [function
  reference](https://adrientaudiere.github.io/bootpq/reference/index.md)
  documents every argument, including how to control the rarefaction
  depth (`sample.size`) and the quantiles reported (`quantile_prob`).
- `bootpq` works on any `phyloseq` object, so the same workflow applies
  to your own data.
