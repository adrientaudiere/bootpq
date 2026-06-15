# Hill-number group tests across rarefaction permutations for phyloseq objects

################################################################################
#' Test multiple times effect of factor on Hill diversity
#'   with different rarefaction even depth
#'
#' @description
#' \lifecycle{experimental}
#'
#' This reduce the risk of a random drawing of an exceptional situation of an
#' unique rarefaction. For each of `nperm` rarefactions to even depth (using
#' [MiscMetabar::rarefy_even_depth_pq()]), the effect of `fact` on each Hill
#' number is tested with [ggstatsplot::ggbetweenstats()].
#'
#' @param physeq (required) A \code{\link[phyloseq]{phyloseq-class}} object
#'   obtained using the `phyloseq` package.
#' @param fact (required) Name of the factor in `physeq@sam_data` used to plot
#'    different lines.
#' @param q (a vector of integer) The list of q values to compute
#'   the hill number H^q. If Null, no hill number are computed. Default value
#'   compute the Hill number 0 (Species richness), the Hill number 1
#'   (exponential of Shannon Index) and the Hill number 2 (inverse of Simpson
#'   Index). Hill numbers are more appropriate in DNA metabarcoding studies
#'   when `q > 0` (Alberdi & Gilbert, 2019; Calderón-Sanou et al., 2019).
#' @param nperm (int) The number of permutations to perform.
#' @param sample.size (int) A single integer value equal to the number of
#'   reads being simulated, also known as the depth. See
#'   [phyloseq::rarefy_even_depth()] and [MiscMetabar::rarefy_even_depth_pq()].
#' @param verbose (logical). If TRUE, print additional information.
#' @param progress_bar (logical, default TRUE) Do we print progress during
#'   the calculation?
#' @param p_val_signif (float, `[0:1]`) The minimum value of p-value to count a
#'   test as significant in the `prop_signif` result.
#' @param type A character specifying the type of statistical approach
#'   (See [ggstatsplot::ggbetweenstats()] for more details):
#'
#'   - "parametric"
#'   - "nonparametric"
#'   - "robust"
#'   - "bayes"
#'
#' @param ... Additional arguments passed on to [ggstatsplot::ggbetweenstats()]
#'   function.
#' @seealso [ggstatsplot::ggbetweenstats()], [MiscMetabar::hill_pq()],
#'   [adonis_rarperm_pq()], [var_par_rarperm_pq()]
#' @return A list of 6 components :
#'
#' - method
#' - expressions
#' - plots
#' - pvals
#' - prop_signif
#' - statistics
#'
#' @export
#' @author Adrien Taudière
#' @references
#' Alberdi, A., & Gilbert, M. T. P. (2019). A guide to the application of
#'   Hill numbers to DNA-based diversity analyses. *Molecular Ecology Resources*.
#'   \doi{10.1111/1755-0998.13014}
#'
#' Calderón-Sanou, I., Münkemüller, T., Boyer, F., Zinger, L., & Thuiller, W.
#'   (2019). From environmental DNA sequences to ecological conclusions: How
#'   strong is the influence of methodological choices? *Journal of Biogeography*,
#'   47. \doi{10.1111/jbi.13681}
#' @examples
#' \dontrun{
#' if (requireNamespace("ggstatsplot")) {
#'   hill_test_rarperm_pq(data_fungi, "Time", nperm = 3)
#'   res <- hill_test_rarperm_pq(data_fungi, "Height",
#'     nperm = 3,
#'     p_val_signif = 0.9
#'   )
#'   patchwork::wrap_plots(res$plots[[1]])
#'   res$plots[[1]][[1]] + res$plots[[2]][[1]] + res$plots[[3]][[1]]
#'   res$prop_signif
#'   res_para <- hill_test_rarperm_pq(data_fungi, "Height",
#'     nperm = 3,
#'     type = "parametric"
#'   )
#'   res_para$plots[[1]][[1]] + res_para$plots[[2]][[1]] + res_para$plots[[3]][[1]]
#'   res_para$pvals
#'   res_para$method
#'   res_para$expressions[[1]]
#' }
#' }
hill_test_rarperm_pq <- function(
  physeq,
  fact,
  q = c(0, 1, 2),
  nperm = 99,
  sample.size = min(phyloseq::sample_sums(physeq)),
  verbose = FALSE,
  progress_bar = TRUE,
  p_val_signif = 0.05,
  type = "nonparametric",
  ...
) {
  if (!requireNamespace("ggstatsplot", quietly = TRUE)) {
    cli::cli_abort(
      "Package {.pkg ggstatsplot} is required for {.fn hill_test_rarperm_pq}. Please install it."
    )
  }
  MiscMetabar::verify_pq(physeq)

  if (nlevels(as.factor(physeq@sam_data[[fact]])) < 2) {
    stop(
      "The factor '",
      fact,
      "' must have at least two levels for ",
      "hill_test_rarperm_pq (statistical tests require at least 2 groups)."
    )
  }
  res_perm <- vector("list", nperm) # pre-allocated for performance
  p_perm <- vector("list", nperm) # pre-allocated for performance
  if (progress_bar) {
    pb <- utils::txtProgressBar(
      min = 0,
      max = nperm * length(q),
      style = 3,
      width = 50,
      char = "="
    )
  }
  if (!exists(".Random.seed", envir = .GlobalEnv)) {
    set.seed(NULL)
  }
  for (i in 1:nperm) {
    if (verbose) {
      psm <-
        MiscMetabar::psmelt_samples_pq(
          physeq = MiscMetabar::rarefy_even_depth_pq(
            physeq,
            rngseed = i,
            sample_size = sample.size
          ),
          q = q
        )
    } else {
      psm <-
        suppressMessages(MiscMetabar::psmelt_samples_pq(
          physeq = MiscMetabar::rarefy_even_depth_pq(
            physeq,
            rngseed = i,
            sample_size = sample.size
          ),
          q = q
        ))
    }
    p_perm[[i]] <- vector("list", length(q))
    res_perm[[i]] <- vector("list", length(q))
    for (j in seq_along(q)) {
      p_perm[[i]][[j]] <-
        ggstatsplot::ggbetweenstats(
          psm,
          !!fact,
          !!paste0("Hill_", q[[j]]),
          type = type,
          ...
        )
      res_perm[[i]][[j]] <-
        ggstatsplot::extract_stats(p_perm[[i]][[j]])
    }
    if (progress_bar) {
      utils::setTxtProgressBar(pb, i * length(q))
    }
  }

  method <- res_perm[[1]][[1]]$subtitle_data[, c(
    "method",
    "effectsize",
    "conf.method"
  )]

  expressions <- sapply(res_perm, function(x) {
    sapply(x, function(xx) {
      xx$subtitle_data$expression
    })
  })
  rownames(expressions) <- paste0("Hill_", q)
  colnames(expressions) <- paste0("ngseed", 1:nperm)

  statistics <- sapply(res_perm, function(x) {
    sapply(x, function(xx) {
      xx$subtitle_data$statistic
    })
  })
  rownames(statistics) <- paste0("Hill_", q)
  colnames(statistics) <- paste0("ngseed", 1:nperm)

  pvals <- sapply(res_perm, function(x) {
    sapply(x, function(xx) {
      xx$subtitle_data$p.value
    })
  })
  rownames(pvals) <- paste0("Hill_", q)
  colnames(pvals) <- paste0("ngseed_", 1:nperm)

  prop_signif <- rowSums(pvals < p_val_signif) / ncol(pvals)
  names(prop_signif) <- paste0("Hill_", q)
  res <-
    list(
      "method" = method,
      "expressions" = expressions,
      "plots" = p_perm,
      "pvals" = pvals,
      "prop_signif" = prop_signif,
      "statistics" = statistics
    )
  return(res)
}
################################################################################
