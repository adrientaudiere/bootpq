# bootpq

`bootpq` provides bootstrapping tools for
[`phyloseq`](https://joey711.github.io/phyloseq/) objects based on
**repeated rarefaction to even sequencing depth**. Instead of drawing
conclusions from a single, possibly atypical, rarefaction, statistical
analyses are computed across many rarefaction permutations — each with a
different random seed — and summarized (mean and quantiles).

`bootpq` is part of the [pqverse](https://github.com/adrientaudiere)
ecosystem and builds on top of
[`MiscMetabar`](https://adrientaudiere.github.io/MiscMetabar/).

## Installation

You can install the development version of `bootpq` from GitHub with:

``` r

# install.packages("remotes")
remotes::install_github("adrientaudiere/bootpq")
```

## Functions

| Function | Description |
|----|----|
| [`adonis_rarperm_pq()`](https://adrientaudiere.github.io/bootpq/reference/adonis_rarperm_pq.md) | PERMANOVA across rarefaction permutations |
| [`hill_test_rarperm_pq()`](https://adrientaudiere.github.io/bootpq/reference/hill_test_rarperm_pq.md) | Hill-number group tests across rarefaction permutations |
| [`var_par_rarperm_pq()`](https://adrientaudiere.github.io/bootpq/reference/var_par_rarperm_pq.md) | Variation partitioning across rarefaction permutations |

## Example

``` r

library(bootpq)
library(MiscMetabar)
data(data_fungi_mini)

data_fungi_woNA <- subset_samples(
  data_fungi_mini,
  !is.na(Time) & !is.na(Height)
)

# PERMANOVA summarized over 99 rarefaction permutations
res <- adonis_rarperm_pq(data_fungi_woNA, "Time*Height", na_remove = TRUE)
res$mean
```
