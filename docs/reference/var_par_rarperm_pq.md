# Partition the Variation of a phyloseq object with rarefaction permutations

**\[experimental\]**

This is an extension of the function
[`MiscMetabar::var_par_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/var_par_pq.html).
The main addition is the computation of `nperm` permutations with
rarefaction even depth by sample. Adjusted R squared values are averaged
across permutations and quantiles are reported.

## Usage

``` r
var_par_rarperm_pq(
  physeq,
  list_component,
  dist_method = "bray",
  nperm = 99,
  quantile_prob = 0.975,
  dbrda_computation = FALSE,
  dbrda_signif_pval = 0.05,
  sample.size = min(phyloseq::sample_sums(physeq)),
  verbose = FALSE,
  progress_bar = TRUE
)
```

## Arguments

- physeq:

  (required) A
  [`phyloseq-class`](https://rdrr.io/pkg/phyloseq/man/phyloseq-class.html)
  object obtained using the `phyloseq` package.

- list_component:

  (required) A named list of 2, 3 or four vectors with names from the
  `@sam_data` slot.

- dist_method:

  (default "bray") the distance used. See
  [`phyloseq::distance()`](https://rdrr.io/pkg/phyloseq/man/distance.html)
  for all available distances or run
  [`phyloseq::distanceMethodList()`](https://rdrr.io/pkg/phyloseq/man/distanceMethodList.html).
  For aitchison and robust.aitchison distance,
  [`vegan::vegdist()`](https://vegandevs.github.io/vegan/reference/vegdist.html)
  function is directly used.

- nperm:

  (int) The number of permutations to perform.

- quantile_prob:

  (float, `[0:1]`) the value to compute the quantile. Minimum quantile
  is compute using `1 - quantile_prob`.

- dbrda_computation:

  (logical) Are dbrda computations run for each individual component
  (each name of the list component) ?

- dbrda_signif_pval:

  (float, `[0:1]`) The value under which the dbrda is considered
  significant.

- sample.size:

  (int) A single integer value equal to the number of reads being
  simulated, also known as the depth. See
  [`phyloseq::rarefy_even_depth()`](https://rdrr.io/pkg/phyloseq/man/rarefy_even_depth.html)
  and
  [`MiscMetabar::rarefy_even_depth_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/rarefy_even_depth_pq.html).

- verbose:

  (logical). If TRUE, print additional information.

- progress_bar:

  (logical, default TRUE) Do we print progress during the calculation?

## Value

A list of class varpart with additional information in the
`$part$indfract` part. Adj.R.square is the mean across permutation.
Adj.R.squared_quantil_min and Adj.R.squared_quantil_max represent the
quantile values of adjusted R squared.

## Details

This function is mainly a wrapper of the work of others. Please make a
reference to
[`vegan::varpart()`](https://vegandevs.github.io/vegan/reference/varpart.html)
if you use this function.

## See also

[`MiscMetabar::var_par_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/var_par_pq.html),
[`vegan::varpart()`](https://vegandevs.github.io/vegan/reference/varpart.html),
[`MiscMetabar::plot_var_part_pq()`](https://adrientaudiere.github.io/MiscMetabar/reference/plot_var_part_pq.html),
[`adonis_rarperm_pq()`](https://adrientaudiere.github.io/bootpq/reference/adonis_rarperm_pq.md),
[`hill_test_rarperm_pq()`](https://adrientaudiere.github.io/bootpq/reference/hill_test_rarperm_pq.md)

## Author

Adrien Taudière

## Examples

``` r
# \donttest{
if (requireNamespace("vegan")) {
  data_fungi_woNA <- subset_samples(
    data_fungi_mini,
    !is.na(Time) & !is.na(Height)
  )
  res_var_2 <- var_par_rarperm_pq(
    data_fungi_woNA,
    list_component = list(
      "Time" = c("Time"),
      "Size" = c("Height", "Diameter")
    ),
    nperm = 2,
    dbrda_computation = TRUE
  )
}
#> Taxa are now in columns.
#>   |                                                          |                                                  |   0%  |                                                          |=========================                         |  50%  |                                                          |==================================================| 100%
# }
if (FALSE) { # \dontrun{
MiscMetabar::plot_var_part_pq(res_var_2)
} # }
```
