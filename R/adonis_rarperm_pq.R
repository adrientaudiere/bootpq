# PERMANOVA across rarefaction permutations for phyloseq objects

################################################################################
#' Permanova (adonis) on permutations of rarefaction even depth
#'
#' @description
#' \lifecycle{experimental}
#'
#' Permanova are computed on a given number of rarefaction with different
#' seed number. This reduce the risk of a random drawing of an exceptional
#' situation of an unique rarefaction. Each permutation rarefies `physeq` to
#' even depth (using [MiscMetabar::rarefy_even_depth_pq()]) and runs
#' [MiscMetabar::adonis_pq()]; results are summarized across permutations.
#'
#' @param physeq (required) A \code{\link[phyloseq]{phyloseq-class}} object
#'   obtained using the `phyloseq` package.
#' @param formula (required) the right part of a formula for
#'   [vegan::adonis2()]. Variables must be present in the `physeq@sam_data`
#'   slot.
#' @param dist_method (default "bray") the distance used. See
#'   [phyloseq::distance()] for all available distances or run
#'   [phyloseq::distanceMethodList()]. For aitchison and robust.aitchison
#'   distance, [vegan::vegdist()] function is directly used.
#' @param merge_sample_by a vector to determine which samples to merge using
#'   the [MiscMetabar::merge_samples2()] function. Need to be in
#'   `physeq@sam_data`.
#' @param na_remove (logical, default FALSE) If set to TRUE, remove samples with
#'   NA in the variables set in formula.
#' @param rarefy_nb_seqs (logical, default FALSE) Rarefy each sample
#'   (before merging if merge_sample_by is set) using
#'   [MiscMetabar::rarefy_even_depth_pq()] inside [MiscMetabar::adonis_pq()].
#' @param verbose (logical, default TRUE) If TRUE, prompt some messages.
#' @param nperm (int, default 99) The number of permutations to perform.
#' @param progress_bar (logical, default TRUE) Do we print progress during
#'   the calculation?
#' @param quantile_prob (float, `[0:1]`) the value to compute the quantile.
#'   Minimum quantile is computed using `1 - quantile_prob`.
#' @param sample.size (int) A single integer value equal to the number of
#'   reads being simulated, also known as the depth. See
#'   [phyloseq::rarefy_even_depth()] and [MiscMetabar::rarefy_even_depth_pq()].
#' @param ... Other params to be passed on to [MiscMetabar::adonis_pq()]
#'   function.
#'
#' @return A list of three dataframe representing the mean, the minimum quantile
#'  and the maximum quantile value for adonis results. See
#'  [MiscMetabar::adonis_pq()].
#' @export
#' @author Adrien Taudière
#' @seealso [MiscMetabar::adonis_pq()], [hill_test_rarperm_pq()],
#'   [var_par_rarperm_pq()]
#' @examples
#' \donttest{
#' if (requireNamespace("vegan")) {
#'   data_fungi_woNA <-
#'     subset_samples(data_fungi_mini, !is.na(Time) & !is.na(Height))
#'   adonis_rarperm_pq(data_fungi_woNA, "Time*Height", na_remove = TRUE, nperm = 3)
#' }
#' }
adonis_rarperm_pq <- function(
  physeq,
  formula,
  dist_method = "bray",
  merge_sample_by = NULL,
  na_remove = FALSE,
  rarefy_nb_seqs = FALSE,
  verbose = TRUE,
  nperm = 99,
  progress_bar = TRUE,
  quantile_prob = 0.975,
  sample.size = min(phyloseq::sample_sums(physeq)),
  ...
) {
  # Ensure .Random.seed exists before the rarefaction below saves/restores it.
  # In a fresh R session .Random.seed is absent until the first random operation.
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    sample.int(1L)
  }

  res_perm <- vector("list", nperm)
  if (progress_bar) {
    pb <- utils::txtProgressBar(
      min = 0,
      max = nperm,
      style = 3,
      width = 50,
      char = "="
    )
  }
  for (i in 1:nperm) {
    res_perm[[i]] <-
      MiscMetabar::adonis_pq(
        MiscMetabar::rarefy_even_depth_pq(
          physeq,
          rngseed = i,
          sample_size = sample.size
        ),
        formula,
        dist_method = dist_method,
        merge_sample_by = merge_sample_by,
        na_remove = na_remove,
        correction_for_sample_size = FALSE,
        rarefy_nb_seqs = rarefy_nb_seqs,
        verbose = verbose,
        sample.size = sample.size,
        ...
      )
    if (progress_bar) {
      utils::setTxtProgressBar(pb, i)
    }
  }

  res_dim <- dim(as.data.frame(res_perm[[1]]))
  res_arr <- array(unlist(res_perm), c(res_dim, nperm))
  res_rownames <- rownames(res_perm[[1]])
  res_colnames <- colnames(res_perm[[1]])

  res_adonis <- list(
    "mean" = apply(res_arr, c(1, 2), mean),
    "quantile_min" = apply(
      res_arr,
      c(1, 2),
      stats::quantile,
      na.rm = TRUE,
      probs = 1 - quantile_prob
    ),
    "quantile_max" = apply(
      res_arr,
      c(1, 2),
      stats::quantile,
      na.rm = TRUE,
      probs = quantile_prob
    )
  )
  for (slot in c("mean", "quantile_min", "quantile_max")) {
    rownames(res_adonis[[slot]]) <- res_rownames
    colnames(res_adonis[[slot]]) <- res_colnames
  }

  return(res_adonis)
}
################################################################################
