skip_on_cran()
library(MiscMetabar)
data(data_fungi_mini)

test_that("adonis_rarperm_pq summarises adonis across permutations", {
  skip_if_not_installed("vegan")
  ps <- subset_samples(data_fungi_mini, !is.na(Time) & !is.na(Height))
  ps <- prune_samples(sample_names(ps)[1:20], ps)
  ps <- clean_pq(ps, silent = TRUE)
  res <- adonis_rarperm_pq(
    ps,
    "Time*Height",
    na_remove = TRUE,
    nperm = 3,
    verbose = FALSE,
    progress_bar = FALSE
  )
  expect_type(res, "list")
  expect_named(res, c("mean", "quantile_min", "quantile_max"))
  expect_true(all(vapply(res, is.matrix, logical(1))))
})
