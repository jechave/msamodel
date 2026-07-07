# MSA prediction layer — posterior prediction from a fitted model.
# Input = a fit object (class "msa_agq"); each predictor node-weights a forward
# model function (calculate_*, in R/model.R) over the stored quadrature nodes to
# return a posterior mean + credible band. @family prediction.

#' Predicted lrmsd_i profile with posterior credible bands from an AGQ fit
#'
#' Propagates the `(a1, a2)` posterior of an adaptive Gauss-Hermite quadrature fit
#' ([fit_lrmsd_i_msa_agq()]) to the predicted per-site divergence profile, returning
#' a posterior mean and a central credible interval at every residue. For each
#' quadrature node the full-model profile `lrmsd_i_msa` is evaluated with
#' [calculate_lrmsd_i_msa()]; the per-site posterior is then those node
#' values weighted by their (normalized) posterior masses `exp(log_weight)`.
#' Everything -- mean and band -- comes from the same quadrature, making no Gaussian
#' assumption (consistent with the fit's parameter summaries). Deterministic.
#'
#' The band is read off the node masses as weighted quantiles, so its tails are only
#' as fine as the node grid. To sharpen the band (especially at high `level`, e.g.
#' 0.99), refit [fit_lrmsd_i_msa_agq()] with a larger `n_nodes`; on the bundled
#' `znb_profile` the band is already stable at the default.
#'
#' @param object An `"msa_agq"` fit, the output of [fit_lrmsd_i_msa_agq()].
#' @param spm_pp Preprocessed data from [preprocess_spm()] (the same used for the fit).
#' @param level Credible-interval coverage (default 0.95).
#' @return A tibble with one row per residue: `i`, `pdb_site`, the posterior summary
#'   of the full-model profile (`lrmsd_i_msa_mean`, `lrmsd_i_msa_lower`,
#'   `lrmsd_i_msa_upper`), and their mean-centred counterparts (`nlrmsd_i_msa_mean`,
#'   `nlrmsd_i_msa_lower`, `nlrmsd_i_msa_upper`).
#' @seealso [fit_lrmsd_i_msa_agq()] (produces the posterior);
#'   [calculate_lrmsd_i_msa()] (evaluated at each node).
#' @family prediction
#' @examples
#' \dontrun{
#' pp  <- preprocess_spm(znb_spm)
#' agq <- fit_lrmsd_i_msa_agq(pp, znb_profile)
#' head(predict_lrmsd_i_msa_agq(agq, pp))
#' }
#' @export
predict_lrmsd_i_msa_agq <- function(object, spm_pp, level = 0.95) {
  if (!inherits(object, "msa_agq")) {
    stop("object must be an 'msa_agq' fit (from fit_lrmsd_i_msa_agq)")
  }
  if (length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1)")
  }
  nodes <- object$nodes
  w <- exp(nodes$log_weight)
  w <- w / sum(w)
  alpha <- (1 - level) / 2

  # lrmsd_i_msa at every node: [site x node] matrix, sites carried from node 1.
  # site_map recovers pdb_site (the forward rung returns only i); left_join keeps
  # the rung's row order (unique i), so first's rows stay aligned with P's rows.
  first <- calculate_lrmsd_i_msa(spm_pp, nodes$a1[1], nodes$a2[1]) %>%
    dplyr::left_join(spm_pp$site_map, by = "i")
  P <- vapply(seq_len(nrow(nodes)), function(k)
    calculate_lrmsd_i_msa(spm_pp, nodes$a1[k], nodes$a2[k])$lrmsd_i_msa,
    numeric(nrow(first)))

  mean_i  <- as.numeric(P %*% w)
  lower_i <- apply(P, 1L, function(row) weighted_quantile(row, w, alpha))
  upper_i <- apply(P, 1L, function(row) weighted_quantile(row, w, 1 - alpha))

  out <- tibble(
    i = first$i, pdb_site = first$pdb_site,
    lrmsd_i_msa_mean = mean_i,
    lrmsd_i_msa_lower = lower_i,
    lrmsd_i_msa_upper = upper_i
  )
  # Mean-centre by a single shift (the mean profile's mean) so the band stays a band,
  # for comparison with mean-centred observed profiles.
  shift <- mean(out$lrmsd_i_msa_mean)
  out %>%
    mutate(
      nlrmsd_i_msa_mean  = lrmsd_i_msa_mean  - shift,
      nlrmsd_i_msa_lower = lrmsd_i_msa_lower - shift,
      nlrmsd_i_msa_upper = lrmsd_i_msa_upper - shift
    )
}

# Node masses of an AGQ fit, normalized to sum 1. Shared by the banded predictors.
#' @noRd
agq_node_weights <- function(nodes) {
  w <- exp(nodes$log_weight)
  w / sum(w)
}

# Posterior mean + central credible band of one per-site quantity across AGQ nodes.
# `values` is a [row x node] matrix (one column per node); returns a tibble with
# <name>_mean/_lower/_upper columns. The band is weighted quantiles of the node masses.
#' @noRd
agq_band <- function(values, w, alpha, name) {
  cols <- list(
    as.numeric(values %*% w),
    apply(values, 1L, function(r) weighted_quantile(r, w, alpha)),
    apply(values, 1L, function(r) weighted_quantile(r, w, 1 - alpha))
  )
  names(cols) <- paste0(name, c("_mean", "_lower", "_upper"))
  tibble::as_tibble(cols)
}

#' Nested-model divergence profiles with posterior credible bands from an AGQ fit
#'
#' Propagates the `(a1, a2)` posterior of an adaptive Gauss-Hermite quadrature fit
#' ([fit_lrmsd_i_msa_agq()]) to all four nested model variants (MM, MS, MA, MSA),
#' returning a posterior mean and a central credible band for each, at every residue.
#' It is [predict_lrmsd_i_msa_agq()] widened from the single full model to the four
#' variants: at each quadrature node the four profiles are evaluated with
#' [calculate_lrmsd_i_nested_models()], then each is summarized across nodes by their
#' (normalized) posterior masses `exp(log_weight)`. All-quadrature, no Gaussian
#' assumption, deterministic.
#'
#' The four variants turn the two selection strengths on or off: MM `(0, 0)`,
#' MS `(a1, 0)`, MA `(0, a2)`, MSA `(a1, a2)`. Bands are read off the node masses as
#' weighted quantiles, so their tails are only as fine as the node grid; to sharpen
#' them, refit with a larger `n_nodes`.
#'
#' @param object An `"msa_agq"` fit, the output of [fit_lrmsd_i_msa_agq()].
#' @param spm_pp Preprocessed data from [preprocess_spm()] (the same used for the fit).
#' @param level Credible-interval coverage (default 0.95).
#' @return A tibble with one row per residue: `i`, `pdb_site`, and a mean/lower/upper
#'   triple for each of the four variants -- `lrmsd_i_mm_{mean,lower,upper}`,
#'   `lrmsd_i_ms_{mean,lower,upper}`, `lrmsd_i_ma_{mean,lower,upper}`,
#'   `lrmsd_i_msa_{mean,lower,upper}` (12 band columns).
#' @seealso [predict_lrmsd_i_msa_agq()] (the single full-model profile);
#'   [calculate_lrmsd_i_nested_models()] (evaluated at each node);
#'   [predict_decomposition_i_msa_agq()] (the phi decomposition, banded).
#' @family prediction
#' @examples
#' \dontrun{
#' pp  <- preprocess_spm(znb_spm)
#' agq <- fit_lrmsd_i_msa_agq(pp, znb_profile)
#' head(predict_lrmsd_i_nested_models_agq(agq, pp))
#' }
#' @export
predict_lrmsd_i_nested_models_agq <- function(object, spm_pp, level = 0.95) {
  if (!inherits(object, "msa_agq")) {
    stop("object must be an 'msa_agq' fit (from fit_lrmsd_i_msa_agq)")
  }
  if (length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1)")
  }
  nodes <- object$nodes
  w     <- agq_node_weights(nodes)
  alpha <- (1 - level) / 2

  # The forward fn carries i + pdb_site itself; take keys from node 1.
  first    <- calculate_lrmsd_i_nested_models(spm_pp, nodes$a1[1], nodes$a2[1])
  variants <- c("lrmsd_i_mm", "lrmsd_i_ms", "lrmsd_i_ma", "lrmsd_i_msa")

  # [site x node] matrix of each variant across nodes, then band it.
  bands <- lapply(variants, function(v) {
    P <- vapply(seq_len(nrow(nodes)), function(k)
      calculate_lrmsd_i_nested_models(spm_pp, nodes$a1[k], nodes$a2[k])[[v]],
      numeric(nrow(first)))
    agq_band(P, w, alpha, v)
  })

  dplyr::bind_cols(first[c("i", "pdb_site")], bands)
}

#' Divergence-decomposition contributions with credible bands from an AGQ fit
#'
#' Propagates the `(a1, a2)` posterior of an adaptive Gauss-Hermite quadrature fit
#' ([fit_lrmsd_i_msa_agq()]) to the site-wise divergence decomposition, returning a
#' posterior mean and a central credible band for each of the three contributions
#' (`phi_mut`, `phi_stab`, `phi_act`), at every residue. At each quadrature node the
#' contributions are evaluated with [calculate_decomposition_i_msa()], then each is
#' summarized across nodes by their (normalized) posterior masses `exp(log_weight)`.
#' All-quadrature, no Gaussian assumption, deterministic -- the decomposition-band
#' path: credible bands on the phi contributions without any sampling.
#'
#' Bands are read off the node masses as weighted quantiles, so their tails are only
#' as fine as the node grid; to sharpen them, refit with a larger `n_nodes`. Note the
#' band on a sum is not the sum of the per-contribution bands (quantiles do not add),
#' but the posterior *means* do add: the three `phi_*_mean` columns sum to the
#' full-model mean profile of [predict_lrmsd_i_msa_agq()].
#'
#' @param object An `"msa_agq"` fit, the output of [fit_lrmsd_i_msa_agq()].
#' @param spm_pp Preprocessed data from [preprocess_spm()] (the same used for the fit).
#' @param level Credible-interval coverage (default 0.95).
#' @return A tibble with one row per residue: `i`, `pdb_site`, and a mean/lower/upper
#'   triple for each contribution -- `phi_mut_{mean,lower,upper}`,
#'   `phi_stab_{mean,lower,upper}`, `phi_act_{mean,lower,upper}` (9 band columns).
#' @seealso [calculate_decomposition_i_msa()] (evaluated at each node);
#'   [predict_lrmsd_i_msa_agq()] (the full-model profile the means sum to);
#'   [predict_lrmsd_i_nested_models_agq()] (the four nested variants, banded).
#' @family prediction
#' @examples
#' \dontrun{
#' pp  <- preprocess_spm(znb_spm)
#' agq <- fit_lrmsd_i_msa_agq(pp, znb_profile)
#' head(predict_decomposition_i_msa_agq(agq, pp))
#' }
#' @export
predict_decomposition_i_msa_agq <- function(object, spm_pp, level = 0.95) {
  if (!inherits(object, "msa_agq")) {
    stop("object must be an 'msa_agq' fit (from fit_lrmsd_i_msa_agq)")
  }
  if (length(level) != 1 || level <= 0 || level >= 1) {
    stop("level must be a single number in (0, 1)")
  }
  nodes <- object$nodes
  w     <- agq_node_weights(nodes)
  alpha <- (1 - level) / 2

  first <- calculate_decomposition_i_msa(spm_pp, nodes$a1[1], nodes$a2[1])
  phis  <- c("phi_mut", "phi_stab", "phi_act")

  bands <- lapply(phis, function(v) {
    P <- vapply(seq_len(nrow(nodes)), function(k)
      calculate_decomposition_i_msa(spm_pp, nodes$a1[k], nodes$a2[k])[[v]],
      numeric(nrow(first)))
    agq_band(P, w, alpha, v)
  })

  dplyr::bind_cols(first[c("i", "pdb_site")], bands)
}
