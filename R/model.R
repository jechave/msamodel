# MSA model forward map (divergence calculators at given (a1, a2))
# Functions to compute the model's predicted structural-divergence profiles

#' Fixation probability of a mutant under the MSA model
#'
#' The MSA model proper: the probability that a single-point mutant fixes under
#' stability selection (strength `a1`) and activity selection (strength `a2`),
#' `p_fix = min(exp(-a1 * ddg), 1) * min(exp(-a2 * ddgact), 1)`. This is a property
#' of a mutant on its own -- it depends only on the mutant's two energy changes and
#' the selection strengths, not on any ensemble -- so it is the elementary quantity
#' an evolutionary-trajectory simulation would evaluate step by step, as well as the
#' primitive the ensemble weights ([weights_jm_spm()]) are built from.
#'
#' Pure and vectorised: `ddg` and `ddgact` may be scalars (one mutant) or
#' equal-length vectors (many mutants), and the result matches their shape. This is
#' *not* normalised -- turning fixation probabilities into averaging weights over a
#' particular ensemble of mutants is a separate, ensemble-specific step.
#'
#' @param ddg Stability free-energy change(s) of the mutant(s) (`ddg_jm`). Scalar or
#'   vector.
#' @param ddgact Activity free-energy change(s) of the mutant(s) (`ddgact_jm`), the
#'   same length as `ddg`.
#' @param a1 Stability selection strength (non-negative). `0` disables stability
#'   selection.
#' @param a2 Activity selection strength (non-negative). `0` disables activity
#'   selection.
#' @return A numeric vector of fixation probabilities, the same length as `ddg`.
#' @seealso [weights_jm_spm()] (normalises these into SPM-ensemble averaging
#'   weights); [calculate_dr2_i_msa()] (the site forward map built on them).
#' @family model
#' @examples
#' \dontrun{
#' # One mutant:
#' pfix_msa(ddg = 1.2, ddgact = 0.4, a1 = 1, a2 = 1)
#' # A whole scan:
#' pp <- preprocess_spm(znb_spm)
#' pfix_msa(pp$energy_data$ddg_jm, pp$energy_data$ddgact_jm, a1 = 1, a2 = 1)
#' }
#' @export
pfix_msa <- function(ddg, ddgact, a1, a2) {
  pstab <- pmin(exp(-a1 * ddg), 1)
  pact  <- pmin(exp(-a2 * ddgact), 1)
  pstab * pact
}

#' Predicted per-site structural divergence at one selection strength
#'
#' Computes the model's predicted per-site structural divergence `dr2_i` (the
#' mean squared C-alpha displacement at each site) for a single pair of selection
#' strengths `(a1, a2)`. Each single-point mutant is weighted by its fixation
#' probability under stability selection (strength `a1`) and activity selection
#' (strength `a2`), and the per-site divergences are averaged over mutants with
#' those weights. This is the elementary prediction the nested-model and fitting
#' functions call repeatedly at different `(a1, a2)`.
#'
#' @param spm_pp Preprocessed single-point-mutation data, the output of
#'   [preprocess_spm()] (the per-mutant energies and the site divergence matrix).
#' @param a1 Stability selection strength (non-negative). `0` disables stability
#'   selection.
#' @param a2 Activity selection strength (non-negative). `0` disables activity
#'   selection.
#' @return A tibble with one row per site: `i` (the site index) and `dr2_i` (the
#'   predicted mean squared displacement at that site).
#' @seealso [calculate_lrmsd_i_nested_models()] which calls this at the four model
#'   variants; [calculate_dr2_n_msa()] for the mode-indexed counterpart.
#' @family model
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' # Full model at a representative selection strength:
#' head(calculate_dr2_i_msa(pp, a1 = 1, a2 = 1))
#' }
#' @export
calculate_dr2_i_msa <- function(spm_pp, a1, a2) {
  dr2_ijm <- spm_pp$dr2_ijm

  # SPM-ensemble weights from the MSA model's fixation probabilities.
  weights_jm <- weights_jm_spm(spm_pp, a1, a2)

  # Weighted average over the mutant (j,m) axis: (dr2_ijm) -> (dr2_i)
  dr2_i <- colSums(dr2_ijm * weights_jm)

  tibble(i = as.integer(colnames(dr2_ijm)), dr2_i = dr2_i)
}

#' Predicted per-site log structural divergence at one selection strength
#'
#' The model's forward prediction in the units the fit works in: the per-site log
#' root-mean-square structural divergence `lrmsd_i = log(sqrt(dr2_i))`, at a single
#' pair of selection strengths `(a1, a2)`. This is the sole owner of the
#' `dr2 -> lrmsd` transform -- the likelihood, the fitters, and the nested-model
#' expansion all obtain their predicted profile through this function rather than
#' re-applying `log(sqrt(.))` by hand.
#'
#' @param spm_pp Preprocessed single-point-mutation data, the output of
#'   [preprocess_spm()].
#' @param a1 Stability selection strength (non-negative). `0` disables stability
#'   selection.
#' @param a2 Activity selection strength (non-negative). `0` disables activity
#'   selection.
#' @return A tibble with one row per site: the site index `i` and the predicted
#'   profile `lrmsd_i_msa`. Callers needing the PDB residue number join
#'   `spm_pp$site_map` themselves (as [calculate_lrmsd_i_nested_models()] does).
#' @seealso [calculate_dr2_i_msa()] (the squared-divergence forward map it wraps);
#'   [calculate_lrmsd_i_nested_models()] (the four-variant analysis expansion).
#' @family model
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' head(calculate_lrmsd_i_msa(pp, a1 = 1, a2 = 1))
#' }
#' @export
calculate_lrmsd_i_msa <- function(spm_pp, a1, a2) {
  dr2 <- calculate_dr2_i_msa(spm_pp, a1, a2)
  tibble(i = dr2$i, lrmsd_i_msa = log(sqrt(dr2$dr2_i)))
}

#' Predicted per-mode structural divergence at one selection strength
#'
#' Mode-indexed counterpart of [calculate_dr2_i_msa()]: instead of a divergence
#' per residue, it returns the predicted divergence `dr2_n` carried by each normal
#' mode of the structure. The mutant fixation weights are exactly those of the
#' site form (they depend only on the mutant, not on whether the output is indexed
#' by site or mode), so the two functions share the same weighting and differ only
#' in what the columns represent. This function predicts only; the package ships no
#' empirical mode profile to fit against.
#'
#' @param spm_pp Preprocessed single-point-mutation data in mode form, the output
#'   of [preprocess_spm_mode()] (the per-mutant energies and the mode divergence
#'   matrix).
#' @param a1 Stability selection strength (non-negative). `0` disables stability
#'   selection.
#' @param a2 Activity selection strength (non-negative). `0` disables activity
#'   selection.
#' @return A tibble with one row per mode: `n` (the mode index) and `dr2_n` (the
#'   predicted mean squared divergence carried by that mode).
#' @seealso [calculate_dr2_i_msa()] for the residue-indexed form;
#'   [calculate_lrmsd_n_nested_models()] which calls this at the four model variants.
#' @family model
#' @examples
#' \dontrun{
#' pp <- preprocess_spm_mode(znb_spm)
#' head(calculate_dr2_n_msa(pp, a1 = 1, a2 = 1))
#' }
#' @export
calculate_dr2_n_msa <- function(spm_pp, a1, a2) {
  dr2_njm <- spm_pp$dr2_njm

  # SPM-ensemble weights from the MSA model's fixation probabilities
  # (axis-agnostic: identical weights as the site form).
  weights_jm <- weights_jm_spm(spm_pp, a1, a2)

  # Weighted average over the mutant (j,m) axis: (dr2_njm) -> (dr2_n)
  dr2_n <- colSums(dr2_njm * weights_jm)

  tibble(n = as.integer(colnames(dr2_njm)), dr2_n = dr2_n)
}

#' Predicted per-mode log structural divergence at one selection strength
#'
#' Mode-indexed counterpart of [calculate_lrmsd_i_msa()]: the model's forward
#' prediction in the units the fit works in, the per-mode log root-mean-square
#' structural divergence `lrmsd_n = log(sqrt(dr2_n))`, at a single pair of selection
#' strengths `(a1, a2)`. This is the sole owner of the mode `dr2 -> lrmsd` transform
#' -- the mode likelihood, the mode fitter, and the mode nested-model expansion all
#' obtain their predicted profile through this function rather than re-applying
#' `log(sqrt(.))` by hand.
#'
#' @param spm_pp Preprocessed single-point-mutation data in mode form, the output
#'   of [preprocess_spm_mode()].
#' @param a1 Stability selection strength (non-negative). `0` disables stability
#'   selection.
#' @param a2 Activity selection strength (non-negative). `0` disables activity
#'   selection.
#' @return A tibble with one row per mode: the mode index `n` and the predicted
#'   profile `lrmsd_n_msa`. Modes are not residue-anchored, so there is no
#'   `pdb_site` column (unlike the site form).
#' @seealso [calculate_dr2_n_msa()] (the squared-divergence forward map it wraps);
#'   [calculate_lrmsd_i_msa()] (the residue-indexed counterpart);
#'   [calculate_lrmsd_n_nested_models()] (the four-variant analysis expansion).
#' @family model
#' @examples
#' \dontrun{
#' pp <- preprocess_spm_mode(znb_spm)
#' head(calculate_lrmsd_n_msa(pp, a1 = 1, a2 = 1))
#' }
#' @export
calculate_lrmsd_n_msa <- function(spm_pp, a1, a2) {
  dr2 <- calculate_dr2_n_msa(spm_pp, a1, a2)
  tibble(n = dr2$n, lrmsd_n_msa = log(sqrt(dr2$dr2_n)))
}

#' Per-site divergence profiles under all four model variants (MM, MS, MA, MSA)
#'
#' The MSA model has two selection pressures, stability and activity, each of
#' which can be switched on or off. This gives four nested variants, from no
#' selection to both, which this function evaluates in one call. For each variant
#' it returns the per-site `lrmsd` profile -- the log root-mean-square structural
#' divergence, `lrmsd = log(sqrt(dr2_i))` -- so the four profiles can be compared
#' or decomposed.
#'
#' The four variants are obtained by turning the two selection strengths on or off:
#'
#' | variant | (a1, a2) | meaning |
#' |---------|----------|---------|
#' | MM  | (0, 0)   | mutation only, no selection |
#' | MS  | (a1, 0)  | mutation + stability selection |
#' | MA  | (0, a2)  | mutation + activity selection |
#' | MSA | (a1, a2) | full model, both selections |
#'
#' Pass the four returned columns to `calculate_msa_decomposition()` to split the
#' divergence into its mutation, stability, and activity contributions.
#'
#' @param spm_pp Preprocessed single-point-mutation data, the output of
#'   [preprocess_spm()].
#' @param a1 Stability selection strength used for the MS and MSA variants
#'   (non-negative).
#' @param a2 Activity selection strength used for the MA and MSA variants
#'   (non-negative).
#' @return A tibble with one row per site: the site index `i`, the PDB residue
#'   number `pdb_site`, and the four divergence profiles `lrmsd_i_mm`,
#'   `lrmsd_i_ms`, `lrmsd_i_ma`, `lrmsd_i_msa`.
#' @seealso [calculate_dr2_i_msa()] (the single-variant prediction it calls);
#'   `calculate_msa_decomposition()` (splits the four columns into contributions);
#'   [calculate_lrmsd_n_nested_models()] (the mode-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' head(calculate_lrmsd_i_nested_models(pp, a1 = 1, a2 = 1))
#' }
#' @export
calculate_lrmsd_i_nested_models <- function(spm_pp, a1, a2) {
  lrmsd <- function(p1, p2) calculate_lrmsd_i_msa(spm_pp, p1, p2)$lrmsd_i_msa

  tibble(
    i           = as.integer(colnames(spm_pp$dr2_ijm)),
    lrmsd_i_mm  = lrmsd(0,  0),
    lrmsd_i_ms  = lrmsd(a1, 0),
    lrmsd_i_ma  = lrmsd(0,  a2),
    lrmsd_i_msa = lrmsd(a1, a2)
  ) %>%
    left_join(spm_pp$site_map, by = "i") %>%
    select(i, pdb_site, everything())
}

#' Per-site divergence decomposition at one (a1, a2)
#'
#' Forward decomposition on the site axis: at a single pair of selection strengths
#' `(a1, a2)` it splits the predicted structural-divergence profile into its
#' mutation, stability, and activity contributions (`phi_mut`, `phi_stab`,
#' `phi_act`). It packages the two-step recipe -- evaluate the four nested model
#' variants with [calculate_lrmsd_i_nested_models()], then apply the sequential
#' split `calculate_msa_decomposition()` -- so a caller wanting phi at a point need
#' not run and thread the four variant columns by hand. The three contributions sum
#' exactly to the full-model profile `lrmsd_i_msa`.
#'
#' This is the single-point (forward) form. Propagating a fitted posterior over
#' `(a1, a2)` to phi with credible bands is done downstream by evaluating this
#' function across the posterior's nodes or draws.
#'
#' @param spm_pp Preprocessed single-point-mutation data, the output of
#'   [preprocess_spm()].
#' @param a1 Stability selection strength (non-negative).
#' @param a2 Activity selection strength (non-negative).
#' @return A tibble with one row per site: the site index `i`, the PDB residue
#'   number `pdb_site`, and the three contributions `phi_mut`, `phi_stab`,
#'   `phi_act`.
#' @seealso [calculate_lrmsd_i_nested_models()] (the four variants it splits);
#'   `calculate_msa_decomposition()` (the pure sequential split);
#'   [calculate_decomposition_n_msa()] (the mode-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' head(calculate_decomposition_i_msa(pp, a1 = 1, a2 = 1))
#' }
#' @export
calculate_decomposition_i_msa <- function(spm_pp, a1, a2) {
  nm  <- calculate_lrmsd_i_nested_models(spm_pp, a1, a2)
  phi <- calculate_msa_decomposition(
    nm$lrmsd_i_mm, nm$lrmsd_i_ms, nm$lrmsd_i_ma, nm$lrmsd_i_msa
  )
  dplyr::bind_cols(nm[c("i", "pdb_site")], tibble::as_tibble(phi))
}

#' Per-mode divergence profiles under all four model variants (MM, MS, MA, MSA)
#'
#' Mode-indexed counterpart of [calculate_lrmsd_i_nested_models()]: it returns the
#' same four nested model variants -- from no selection (MM) to both stability and
#' activity selection (MSA) -- but as a divergence profile over the structure's
#' normal modes rather than its residues. For each variant it reports the per-mode
#' `lrmsd`, the log root-mean-square divergence `lrmsd = log(sqrt(dr2_n))`.
#'
#' The four variants are obtained by turning the two selection strengths on or off:
#'
#' | variant | (a1, a2) | meaning |
#' |---------|----------|---------|
#' | MM  | (0, 0)   | mutation only, no selection |
#' | MS  | (a1, 0)  | mutation + stability selection |
#' | MA  | (0, a2)  | mutation + activity selection |
#' | MSA | (a1, a2) | full model, both selections |
#'
#' Modes are not anchored to residues, so the output has no `pdb_site` column
#' (unlike the site form). Pass the four returned columns to
#' `calculate_msa_decomposition()` to split the divergence into its contributions.
#'
#' @param spm_pp Preprocessed single-point-mutation data in mode form, the output
#'   of [preprocess_spm_mode()].
#' @param a1 Stability selection strength used for the MS and MSA variants
#'   (non-negative).
#' @param a2 Activity selection strength used for the MA and MSA variants
#'   (non-negative).
#' @return A tibble with one row per mode: the mode index `n` and the four
#'   divergence profiles `lrmsd_n_mm`, `lrmsd_n_ms`, `lrmsd_n_ma`, `lrmsd_n_msa`.
#' @seealso [calculate_dr2_n_msa()] (the single-variant prediction it calls);
#'   [calculate_lrmsd_i_nested_models()] (the residue-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' pp <- preprocess_spm_mode(znb_spm)
#' head(calculate_lrmsd_n_nested_models(pp, a1 = 1, a2 = 1))
#' }
#' @export
calculate_lrmsd_n_nested_models <- function(spm_pp, a1, a2) {
  lrmsd <- function(p1, p2) calculate_lrmsd_n_msa(spm_pp, p1, p2)$lrmsd_n_msa

  tibble(
    n           = as.integer(colnames(spm_pp$dr2_njm)),
    lrmsd_n_mm  = lrmsd(0,  0),
    lrmsd_n_ms  = lrmsd(a1, 0),
    lrmsd_n_ma  = lrmsd(0,  a2),
    lrmsd_n_msa = lrmsd(a1, a2)
  )
}

#' Per-mode divergence decomposition at one (a1, a2)
#'
#' Mode-indexed counterpart of [calculate_decomposition_i_msa()]: at a single pair
#' of selection strengths `(a1, a2)` it splits the predicted mode-divergence profile
#' into its mutation, stability, and activity contributions (`phi_mut`, `phi_stab`,
#' `phi_act`). It packages the two-step recipe -- evaluate the four nested model
#' variants with [calculate_lrmsd_n_nested_models()], then apply the sequential split
#' `calculate_msa_decomposition()`. The three contributions sum exactly to the
#' full-model profile `lrmsd_n_msa`. Modes are not anchored to residues, so the
#' output has no `pdb_site` column.
#'
#' This is the single-point (forward) form. Propagating a fitted posterior over
#' `(a1, a2)` to phi with credible bands is done downstream by evaluating this
#' function across the posterior's nodes or draws.
#'
#' @param spm_pp Preprocessed single-point-mutation data in mode form, the output
#'   of [preprocess_spm_mode()].
#' @param a1 Stability selection strength (non-negative).
#' @param a2 Activity selection strength (non-negative).
#' @return A tibble with one row per mode: the mode index `n` and the three
#'   contributions `phi_mut`, `phi_stab`, `phi_act`.
#' @seealso [calculate_lrmsd_n_nested_models()] (the four variants it splits);
#'   `calculate_msa_decomposition()` (the pure sequential split);
#'   [calculate_decomposition_i_msa()] (the residue-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' pp <- preprocess_spm_mode(znb_spm)
#' head(calculate_decomposition_n_msa(pp, a1 = 1, a2 = 1))
#' }
#' @export
calculate_decomposition_n_msa <- function(spm_pp, a1, a2) {
  nm  <- calculate_lrmsd_n_nested_models(spm_pp, a1, a2)
  phi <- calculate_msa_decomposition(
    nm$lrmsd_n_mm, nm$lrmsd_n_ms, nm$lrmsd_n_ma, nm$lrmsd_n_msa
  )
  dplyr::bind_cols(nm["n"], tibble::as_tibble(phi))
}

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
#' @family model
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
#' @family model
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
#' All-quadrature, no Gaussian assumption, deterministic. This is the assumption-free
#' replacement for the MCMC decomposition bands.
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
#' @family model
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
