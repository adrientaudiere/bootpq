skip_on_cran()
library(MiscMetabar)
data(data_fungi_mini)

balanced_height_pq <- function(physeq, n = 6) {
  sd <- as(phyloseq::sample_data(physeq), "data.frame")
  sn <- phyloseq::sample_names(physeq)
  sn_low <- sn[which(sd$Height == "Low")][seq_len(n)]
  sn_high <- sn[which(sd$Height == "High")][seq_len(n)]
  clean_pq(prune_samples(c(sn_low, sn_high), physeq), silent = TRUE)
}

test_that("hill_test_rarperm_pq returns the 6-component result", {
  skip_if_not_installed("ggstatsplot")
  ps <- balanced_height_pq(data_fungi_mini, n = 6)
  res <- hill_test_rarperm_pq(
    ps,
    "Height",
    q = c(0, 1),
    nperm = 2,
    progress_bar = FALSE
  )
  expect_named(
    res,
    c("method", "expressions", "plots", "pvals", "prop_signif", "statistics")
  )
  expect_equal(nrow(res$pvals), 2)
  expect_equal(ncol(res$pvals), 2)
})

test_that("hill_test_rarperm_pq errors with a single-level factor", {
  skip_if_not_installed("ggstatsplot")
  sd <- as(phyloseq::sample_data(data_fungi_mini), "data.frame")
  sn <- phyloseq::sample_names(data_fungi_mini)
  ps <- prune_samples(sn[which(sd$Height == "Low")][1:6], data_fungi_mini)
  expect_error(
    hill_test_rarperm_pq(ps, "Height", nperm = 2, progress_bar = FALSE),
    "at least two levels"
  )
})
