skip_on_cran()
library(MiscMetabar)
data(data_fungi_mini)

test_that("var_par_rarperm_pq averages adjusted R squared across permutations", {
  skip_if_not_installed("vegan")
  ps <- subset_samples(data_fungi_mini, !is.na(Time) & !is.na(Height))
  ps <- prune_samples(sample_names(ps)[1:25], ps)
  ps <- clean_pq(ps, silent = TRUE)
  res <- var_par_rarperm_pq(
    ps,
    list_component = list("Time" = c("Time"), "Size" = c("Height", "Diameter")),
    nperm = 2,
    progress_bar = FALSE
  )
  expect_s3_class(res, "varpart")
  expect_equal(res$Xnames, c("Time", "Size"))
  expect_false(is.null(res$part$indfract$Adj.R.squared_quantil_max))
})
