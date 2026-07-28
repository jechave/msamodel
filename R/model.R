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

#' Predicted mean-centred per-site log structural divergence at one selection strength
#'
#' Mean-centred counterpart of [calculate_lrmsd_i_msa()]: the per-site profile centred by
#' its own mean over the **full model support** (all model residues),
#' `nlrmsd_i = lrmsd_i - mean_S(lrmsd_i)`. This is the quantity the site ML fit is on -- the
#' likelihood centres both sides -- so it is the natural forward map for exploring a
#' scenario `(a1, a2)` on the same scale the fit reports, without fitting. It is also the
#' point profile that [predict_nlrmsd_i_msa_ml()] bands: the two agree exactly (the
#' predictor's `_mean` column is this function's `nlrmsd_i_msa`).
#'
#' Centring is over all model residues, deliberately agnostic to which residues a dataset
#' happens to observe. The fit instead centres over the observation-matched residues, so a
#' fitted profile sits at a slightly different level by design; see [predict_nlrmsd_i_msa_ml()].
#'
#' @param spm_pp Preprocessed single-point-mutation data, the output of [preprocess_spm()].
#' @param a1 Stability selection strength (non-negative). `0` disables stability selection.
#' @param a2 Activity selection strength (non-negative). `0` disables activity selection.
#' @return A tibble with one row per site: the site index `i` and the mean-centred profile
#'   `nlrmsd_i_msa`.
#' @seealso [calculate_lrmsd_i_msa()] (the uncentred profile it centres);
#'   [predict_nlrmsd_i_msa_ml()] (the same quantity with delta-method bands from a fit);
#'   [calculate_nlrmsd_n_msa()] (the mode-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' head(calculate_nlrmsd_i_msa(pp, a1 = 1, a2 = 1))
#' }
#' @export
calculate_nlrmsd_i_msa <- function(spm_pp, a1, a2) {
  lrmsd <- calculate_lrmsd_i_msa(spm_pp, a1, a2)
  tibble(i = lrmsd$i, nlrmsd_i_msa = lrmsd$lrmsd_i_msa - mean(lrmsd$lrmsd_i_msa))
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

#' Predicted mean-centred per-mode log structural divergence at one selection strength
#'
#' Mode-indexed counterpart of [calculate_nlrmsd_i_msa()]: the per-mode profile centred by
#' its own mean over the full model support (all modes),
#' `nlrmsd_n = lrmsd_n - mean_S(lrmsd_n)`. This is the quantity the mode ML fit is on and the
#' point profile that [predict_nlrmsd_n_msa_ml()] bands (the two agree exactly). Modes are
#' not residue-anchored, so there is no `pdb_site` column. See [calculate_nlrmsd_i_msa()] for
#' the support note (predict/scenario centres over all modes; the fit centres over the
#' observation-matched modes, so the two levels differ by design).
#'
#' @param spm_pp Preprocessed single-point-mutation data in mode form, the output of
#'   [preprocess_spm_mode()].
#' @param a1 Stability selection strength (non-negative). `0` disables stability selection.
#' @param a2 Activity selection strength (non-negative). `0` disables activity selection.
#' @return A tibble with one row per mode: the mode index `n` and the mean-centred profile
#'   `nlrmsd_n_msa`.
#' @seealso [calculate_lrmsd_n_msa()] (the uncentred profile it centres);
#'   [predict_nlrmsd_n_msa_ml()] (the same quantity with delta-method bands from a fit);
#'   [calculate_nlrmsd_i_msa()] (the residue-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' pp <- preprocess_spm_mode(znb_spm)
#' head(calculate_nlrmsd_n_msa(pp, a1 = 1, a2 = 1))
#' }
#' @export
calculate_nlrmsd_n_msa <- function(spm_pp, a1, a2) {
  lrmsd <- calculate_lrmsd_n_msa(spm_pp, a1, a2)
  tibble(n = lrmsd$n, nlrmsd_n_msa = lrmsd$lrmsd_n_msa - mean(lrmsd$lrmsd_n_msa))
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

#' Mean-centred per-site divergence profiles under all four model variants (MM, MS, MA, MSA)
#'
#' Mean-centred counterpart of [calculate_lrmsd_i_nested_models()]: the same four nested
#' variants, but each centred by **its own** mean over the full model support (all model
#' residues), `nlrmsd_i_v = lrmsd_i_v - mean_S(lrmsd_i_v)`. Each variant is centred
#' independently -- MM by the MM mean, MSA by the MSA mean -- matching how the fit and
#' [predict_nlrmsd_i_nested_models_ml()] treat them. This is the natural scenario forward map
#' for comparing the centred variants without fitting.
#'
#' @param spm_pp Preprocessed single-point-mutation data, the output of [preprocess_spm()].
#' @param a1 Stability selection strength used for the MS and MSA variants (non-negative).
#' @param a2 Activity selection strength used for the MA and MSA variants (non-negative).
#' @return A tibble with one row per site: `i`, `pdb_site`, and the four mean-centred profiles
#'   `nlrmsd_i_mm`, `nlrmsd_i_ms`, `nlrmsd_i_ma`, `nlrmsd_i_msa`.
#' @seealso [calculate_lrmsd_i_nested_models()] (the uncentred variants it centres);
#'   [predict_nlrmsd_i_nested_models_ml()] (the same variants with delta-method bands);
#'   [calculate_nlrmsd_n_nested_models()] (the mode-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' head(calculate_nlrmsd_i_nested_models(pp, a1 = 1, a2 = 1))
#' }
#' @export
calculate_nlrmsd_i_nested_models <- function(spm_pp, a1, a2) {
  nm <- calculate_lrmsd_i_nested_models(spm_pp, a1, a2)
  tibble(
    i           = nm$i,
    pdb_site    = nm$pdb_site,
    nlrmsd_i_mm  = nm$lrmsd_i_mm  - mean(nm$lrmsd_i_mm),
    nlrmsd_i_ms  = nm$lrmsd_i_ms  - mean(nm$lrmsd_i_ms),
    nlrmsd_i_ma  = nm$lrmsd_i_ma  - mean(nm$lrmsd_i_ma),
    nlrmsd_i_msa = nm$lrmsd_i_msa - mean(nm$lrmsd_i_msa)
  )
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

#' Mean-centred per-site divergence decomposition at one (a1, a2)
#'
#' Mean-centred counterpart of [calculate_decomposition_i_msa()]: the three contributions of
#' the *centred* per-site profile `nlrmsd_i_msa`, each centred by its own mean over the full
#' model support and emitted as `nphi_mut`, `nphi_stab`, `nphi_act`. The decomposition is
#' defined on the centred quantity (the one the fit is on); the three columns sum exactly to
#' `nlrmsd_i_msa` (centring is linear). This is the point decomposition that
#' [predict_nlrmsd_i_msa_decomposition_ml()] bands (the two agree exactly).
#'
#' @param spm_pp Preprocessed single-point-mutation data, the output of [preprocess_spm()].
#' @param a1 Stability selection strength (non-negative).
#' @param a2 Activity selection strength (non-negative).
#' @return A tibble with one row per site: `i`, `pdb_site`, and the three mean-centred
#'   contributions `nphi_mut`, `nphi_stab`, `nphi_act`.
#' @seealso [calculate_decomposition_i_msa()] (the uncentred `phi_*` decomposition);
#'   [calculate_nlrmsd_i_msa()] (the centred profile these sum to);
#'   [predict_nlrmsd_i_msa_decomposition_ml()] (the same contributions with delta-method
#'   bands); [calculate_nlrmsd_n_msa_decomposition()] (the mode-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' head(calculate_nlrmsd_i_msa_decomposition(pp, a1 = 1, a2 = 1))
#' }
#' @export
calculate_nlrmsd_i_msa_decomposition <- function(spm_pp, a1, a2) {
  d <- calculate_decomposition_i_msa(spm_pp, a1, a2)
  tibble(
    i         = d$i,
    pdb_site  = d$pdb_site,
    nphi_mut  = d$phi_mut  - mean(d$phi_mut),
    nphi_stab = d$phi_stab - mean(d$phi_stab),
    nphi_act  = d$phi_act  - mean(d$phi_act)
  )
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

#' Mean-centred per-mode divergence profiles under all four model variants (MM, MS, MA, MSA)
#'
#' Mode-indexed counterpart of [calculate_nlrmsd_i_nested_models()]: the four nested variants
#' of the per-mode profile, each centred by its own mean over the full model support (all
#' modes). Modes are not residue-anchored, so there is no `pdb_site` column. Each variant is
#' centred independently, matching [predict_nlrmsd_n_nested_models_ml()].
#'
#' @param spm_pp Preprocessed single-point-mutation data in mode form, the output of
#'   [preprocess_spm_mode()].
#' @param a1 Stability selection strength used for the MS and MSA variants (non-negative).
#' @param a2 Activity selection strength used for the MA and MSA variants (non-negative).
#' @return A tibble with one row per mode: `n` and the four mean-centred profiles
#'   `nlrmsd_n_mm`, `nlrmsd_n_ms`, `nlrmsd_n_ma`, `nlrmsd_n_msa`.
#' @seealso [calculate_lrmsd_n_nested_models()] (the uncentred variants it centres);
#'   [predict_nlrmsd_n_nested_models_ml()] (the same variants with delta-method bands);
#'   [calculate_nlrmsd_i_nested_models()] (the residue-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' pp <- preprocess_spm_mode(znb_spm)
#' head(calculate_nlrmsd_n_nested_models(pp, a1 = 1, a2 = 1))
#' }
#' @export
calculate_nlrmsd_n_nested_models <- function(spm_pp, a1, a2) {
  nm <- calculate_lrmsd_n_nested_models(spm_pp, a1, a2)
  tibble(
    n           = nm$n,
    nlrmsd_n_mm  = nm$lrmsd_n_mm  - mean(nm$lrmsd_n_mm),
    nlrmsd_n_ms  = nm$lrmsd_n_ms  - mean(nm$lrmsd_n_ms),
    nlrmsd_n_ma  = nm$lrmsd_n_ma  - mean(nm$lrmsd_n_ma),
    nlrmsd_n_msa = nm$lrmsd_n_msa - mean(nm$lrmsd_n_msa)
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

#' Mean-centred per-mode divergence decomposition at one (a1, a2)
#'
#' Mode-indexed counterpart of [calculate_nlrmsd_i_msa_decomposition()]: the three
#' contributions of the centred per-mode profile `nlrmsd_n_msa`, each centred by its own mean
#' over the full model support and emitted as `nphi_mut`, `nphi_stab`, `nphi_act`. The three
#' columns sum exactly to `nlrmsd_n_msa`. Modes are not residue-anchored, so there is no
#' `pdb_site` column. This is the point decomposition that
#' [predict_nlrmsd_n_msa_decomposition_ml()] bands.
#'
#' @param spm_pp Preprocessed single-point-mutation data in mode form, the output of
#'   [preprocess_spm_mode()].
#' @param a1 Stability selection strength (non-negative).
#' @param a2 Activity selection strength (non-negative).
#' @return A tibble with one row per mode: `n` and the three mean-centred contributions
#'   `nphi_mut`, `nphi_stab`, `nphi_act`.
#' @seealso [calculate_decomposition_n_msa()] (the uncentred `phi_*` decomposition);
#'   [calculate_nlrmsd_n_msa()] (the centred profile these sum to);
#'   [predict_nlrmsd_n_msa_decomposition_ml()] (the same contributions with delta-method
#'   bands); [calculate_nlrmsd_i_msa_decomposition()] (the residue-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' pp <- preprocess_spm_mode(znb_spm)
#' head(calculate_nlrmsd_n_msa_decomposition(pp, a1 = 1, a2 = 1))
#' }
#' @export
calculate_nlrmsd_n_msa_decomposition <- function(spm_pp, a1, a2) {
  d <- calculate_decomposition_n_msa(spm_pp, a1, a2)
  tibble(
    n         = d$n,
    nphi_mut  = d$phi_mut  - mean(d$phi_mut),
    nphi_stab = d$phi_stab - mean(d$phi_stab),
    nphi_act  = d$phi_act  - mean(d$phi_act)
  )
}
