# Variation partitioning across rarefaction permutations for phyloseq objects

################################################################################
#' Partition the Variation of a phyloseq object with rarefaction permutations
#' @description
#' \lifecycle{experimental}
#'
#'   This is an extension of the function [MiscMetabar::var_par_pq()]. The main
#'   addition is the computation of `nperm` permutations with rarefaction even
#'   depth by sample. Adjusted R squared values are averaged across permutations
#'   and quantiles are reported.
#'
#' @param physeq (required) A \code{\link[phyloseq]{phyloseq-class}} object
#'   obtained using the `phyloseq` package.
#' @param list_component (required) A named list of 2, 3 or four vectors with
#'   names from the `@sam_data` slot.
#' @param dist_method (default "bray") the distance used. See
#'   [phyloseq::distance()] for all available distances or run
#'   [phyloseq::distanceMethodList()].
#'   For aitchison and robust.aitchison distance, [vegan::vegdist()]
#'   function is directly used.
#' @param nperm (int) The number of permutations to perform.
#' @param quantile_prob (float, `[0:1]`) the value to compute the quantile.
#'   Minimum quantile is compute using `1 - quantile_prob`.
#' @param dbrda_computation (logical) Are dbrda computations run for each
#'  individual component (each name of the list component) ?
#' @param dbrda_signif_pval (float, `[0:1]`) The value under which the dbrda is
#'   considered significant.
#' @param sample.size (int) A single integer value equal to the number of
#'   reads being simulated, also known as the depth. See
#'   [phyloseq::rarefy_even_depth()] and [MiscMetabar::rarefy_even_depth_pq()].
#' @param verbose (logical). If TRUE, print additional information.
#' @param progress_bar (logical, default TRUE) Do we print progress during
#'   the calculation?
#'
#' @return A list of class varpart with additional information in the
#'  `$part$indfract` part. Adj.R.square is the mean across permutation.
#'   Adj.R.squared_quantil_min and Adj.R.squared_quantil_max represent
#'   the quantile values of adjusted R squared.
#' @export
#' @seealso [MiscMetabar::var_par_pq()], [vegan::varpart()],
#'   [MiscMetabar::plot_var_part_pq()], [adonis_rarperm_pq()],
#'   [hill_test_rarperm_pq()]
#' @author Adrien Taudière
#' @examples
#' \donttest{
#' if (requireNamespace("vegan")) {
#'   data_fungi_woNA <- subset_samples(
#'     data_fungi_mini,
#'     !is.na(Time) & !is.na(Height)
#'   )
#'   res_var_2 <- var_par_rarperm_pq(
#'     data_fungi_woNA,
#'     list_component = list(
#'       "Time" = c("Time"),
#'       "Size" = c("Height", "Diameter")
#'     ),
#'     nperm = 2,
#'     dbrda_computation = TRUE
#'   )
#' }
#' }
#' \dontrun{
#' MiscMetabar::plot_var_part_pq(res_var_2)
#' }
#' @details
#' This function is mainly a wrapper of the work of others.
#'   Please make a reference to `vegan::varpart()` if you
#'   use this function.
var_par_rarperm_pq <-
  function(
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
  ) {
    physeq <- MiscMetabar::taxa_as_columns(physeq)
    MiscMetabar::verify_pq(physeq)

    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      sample.int(1L)
    }

    if (progress_bar) {
      pb <- utils::txtProgressBar(
        min = 0,
        max = nperm,
        style = 3,
        width = 50,
        char = "="
      )
    }

    res_perm <- vector("list", nperm)
    for (i in 1:nperm) {
      res_perm[[i]] <-
        MiscMetabar::var_par_pq(
          physeq = MiscMetabar::rarefy_even_depth_pq(
            physeq,
            rngseed = i,
            sample_size = sample.size
          ),
          list_component = list_component,
          dist_method = dist_method,
          dbrda_computation = dbrda_computation
        )

      if (progress_bar) {
        utils::setTxtProgressBar(pb, i)
      }
    }
    res_varpart <- MiscMetabar::var_par_pq(
      physeq = physeq,
      list_component = list_component,
      dist_method = dist_method,
      dbrda_computation = dbrda_computation
    )

    if (dbrda_computation) {
      res_varpart$dbrda_result_prop_pval_signif <-
        rowSums(
          sapply(res_perm, function(x) {
            sapply(x$dbrda_result, function(xx) {
              xx$`Pr(>F)`[[1]]
            })
          }) <
            dbrda_signif_pval
        ) /
        nperm
    }

    # Pre-compute sapply results for efficiency
    r_square_matrix <- sapply(res_perm, function(x) x$part$indfract$R.square)
    adj_r_square_matrix <- sapply(res_perm, function(x) {
      x$part$indfract$Adj.R.square
    })

    res_varpart$part$indfract$R.square <- rowMeans(r_square_matrix)
    res_varpart$part$indfract$R.square_quantil_max <-
      apply(r_square_matrix, 1, function(xx) {
        stats::quantile(xx, probs = quantile_prob, na.rm = TRUE)
      })
    res_varpart$part$indfract$R.square_quantil_min <-
      apply(r_square_matrix, 1, function(xx) {
        stats::quantile(xx, probs = 1 - quantile_prob, na.rm = TRUE)
      })

    res_varpart$part$indfract$Adj.R.square <- rowMeans(adj_r_square_matrix)
    res_varpart$part$indfract$Adj.R.squared_quantil_max <-
      apply(adj_r_square_matrix, 1, function(xx) {
        stats::quantile(xx, probs = quantile_prob, na.rm = TRUE)
      })
    res_varpart$part$indfract$Adj.R.squared_quantil_min <-
      apply(adj_r_square_matrix, 1, function(xx) {
        stats::quantile(xx, probs = 1 - quantile_prob, na.rm = TRUE)
      })

    res_varpart$Xnames <- names(list_component)
    return(res_varpart)
  }
################################################################################
