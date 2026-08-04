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
#' spm <- generate_spm_data(znb_wt, seed = 1024)
#' pfix_msa(spm$energy_data$ddg_jm, spm$energy_data$ddgact_jm, a1 = 1, a2 = 1)
#' }
#' @export
pfix_msa <- function(ddg, ddgact, a1, a2) {
  pstab <- pmin(exp(-a1 * ddg), 1)
  pact  <- pmin(exp(-a2 * ddgact), 1)
  pstab * pact
}

# ---- axis-blind forward-map primitives (internal) ---------------------------------
# The site (_i_/dr2_ijm) and mode (_n_/dr2_njm) forward maps are the same math on a
# different response-mutation matrix. These primitives hold that shared math: they take
# a bare [mutant x response] matrix `dr2_mat` (either dr2_ijm or dr2_njm) plus the
# `energy_data` tibble, and return bare vectors (or a named list of bare vectors) in
# matrix-column order -- no index, no tibble, no site_map. Alignment is by column
# POSITION (the dr2 matrices carry no column names by construction); the exported
# calculate_*_i/n skins below attach the axis index (and, site only, the pdb_site label)
# at the boundary.

#' SPM-ensemble averaging weights from an energy table
#'
#' Turns the per-mutant MSA fixation probabilities into the normalised averaging weights
#' `weights_jm = pfix_jm / sum(pfix_jm)`, one per mutant `(j, m)`, summing to one. Reads
#' the per-mutant energies directly from an `energy_data` tibble, so the forward-map
#' primitives obtain their weights without the `spm` object.
#'
#' @param energy_data A tibble carrying per-mutant `ddg_jm` and `ddgact_jm` columns.
#' @param a1 Stability selection strength (non-negative).
#' @param a2 Activity selection strength (non-negative).
#' @seealso [weights_jm_spm()], the exported form taking the whole `spm` object.
#' @return A numeric vector of averaging weights, one per mutant, summing to one.
#' @family model
#' @noRd
weights_jm <- function(energy_data, a1, a2) {
  pfix_jm <- pfix_msa(energy_data$ddg_jm, energy_data$ddgact_jm, a1, a2)
  pfix_jm / sum(pfix_jm)
}

#' Axis-blind per-response structural divergence at one (a1, a2)
#'
#' The shared core of [calculate_dr2_i_msa()] / [calculate_dr2_n_msa()]: weights each
#' mutant by its MSA fixation probability and averages the per-response squared
#' displacements over mutants.
#'
#' @param dr2_mat A `[mutant x response]` divergence matrix (`dr2_ijm` or `dr2_njm`).
#' @param energy_data The per-mutant energy tibble (for the weights).
#' @param a1,a2 Selection strengths.
#' @return A numeric vector of per-response `dr2`, in column order.
#' @family model
#' @noRd
dr2_msa <- function(dr2_mat, energy_data, a1, a2) {
  w <- weights_jm(energy_data, a1, a2)
  colSums(dr2_mat * w)
}

#' Axis-blind per-response log structural divergence at one (a1, a2)
#'
#' The shared core of [calculate_lrmsd_i_msa()] / [calculate_lrmsd_n_msa()]:
#' `lrmsd = log(sqrt(dr2))`. Sole owner of the `dr2 -> lrmsd` transform.
#'
#' @inheritParams dr2_msa
#' @return A numeric vector of per-response `lrmsd`, in column order.
#' @family model
#' @noRd
lrmsd_msa <- function(dr2_mat, energy_data, a1, a2) {
  log(sqrt(dr2_msa(dr2_mat, energy_data, a1, a2)))
}

#' Axis-blind mean-centred per-response log divergence at one (a1, a2)
#'
#' The shared core of [calculate_nlrmsd_i_msa()] / [calculate_nlrmsd_n_msa()]:
#' `nlrmsd = lrmsd - mean(lrmsd)` over the full response support.
#'
#' @inheritParams dr2_msa
#' @return A numeric vector of per-response `nlrmsd`, in column order.
#' @family model
#' @noRd
nlrmsd_msa <- function(dr2_mat, energy_data, a1, a2) {
  lrmsd <- lrmsd_msa(dr2_mat, energy_data, a1, a2)
  lrmsd - mean(lrmsd)
}

#' Axis-blind four nested-model lrmsd profiles at one (a1, a2)
#'
#' The shared core of [calculate_lrmsd_i_nested_models()] /
#' [calculate_lrmsd_n_nested_models()]: the four variants MM `(0,0)`, MS `(a1,0)`,
#' MA `(0,a2)`, MSA `(a1,a2)`, each an `lrmsd` vector.
#'
#' @inheritParams dr2_msa
#' @return A named list `mm`, `ms`, `ma`, `msa`, each an `lrmsd` vector in column order.
#' @family model
#' @noRd
lrmsd_nested_models <- function(dr2_mat, energy_data, a1, a2) {
  lrmsd <- function(p1, p2) lrmsd_msa(dr2_mat, energy_data, p1, p2)
  list(mm  = lrmsd(0,  0),
       ms  = lrmsd(a1, 0),
       ma  = lrmsd(0,  a2),
       msa = lrmsd(a1, a2))
}

#' Axis-blind four mean-centred nested-model nlrmsd profiles at one (a1, a2)
#'
#' The shared core of [calculate_nlrmsd_i_nested_models()] /
#' [calculate_nlrmsd_n_nested_models()]: each nested variant centred by its own mean.
#'
#' @inheritParams dr2_msa
#' @return A named list `mm`, `ms`, `ma`, `msa`, each a centred vector in column order.
#' @family model
#' @noRd
nlrmsd_nested_models <- function(dr2_mat, energy_data, a1, a2) {
  v <- lrmsd_nested_models(dr2_mat, energy_data, a1, a2)
  lapply(v, function(x) x - mean(x))
}

#' Axis-blind divergence decomposition at one (a1, a2)
#'
#' The shared core of [calculate_decomposition_i_msa()] / [calculate_decomposition_n_msa()]:
#' evaluate the four nested variants, then apply the sequential split
#' `decompose_nested()`. Named `lrmsd_msa_decomposition` (not the naive
#' `decomposition_msa`) so it pairs cleanly with `nlrmsd_msa_decomposition`; the exported
#' `calculate_decomposition_*_msa` is a known naming wart kept for now (see plan).
#'
#' @inheritParams dr2_msa
#' @return A list `phi_mut`, `phi_stab`, `phi_act`, each a vector in column order.
#' @family model
#' @noRd
lrmsd_msa_decomposition <- function(dr2_mat, energy_data, a1, a2) {
  v <- lrmsd_nested_models(dr2_mat, energy_data, a1, a2)
  decompose_nested(v$mm, v$ms, v$ma, v$msa)
}

#' Axis-blind mean-centred divergence decomposition at one (a1, a2)
#'
#' The shared core of [calculate_nlrmsd_i_msa_decomposition()] /
#' [calculate_nlrmsd_n_msa_decomposition()]: each phi contribution centred by its own mean.
#'
#' @inheritParams dr2_msa
#' @return A named list `nphi_mut`, `nphi_stab`, `nphi_act`, each a centred vector in column order.
#' @family model
#' @noRd
nlrmsd_msa_decomposition <- function(dr2_mat, energy_data, a1, a2) {
  phi <- lrmsd_msa_decomposition(dr2_mat, energy_data, a1, a2)
  list(nphi_mut  = phi$phi_mut  - mean(phi$phi_mut),
       nphi_stab = phi$phi_stab - mean(phi$phi_stab),
       nphi_act  = phi$phi_act  - mean(phi$phi_act))
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
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
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
#' spm <- generate_spm_data(znb_wt, seed = 1024)
#' # Full model at a representative selection strength:
#' head(calculate_dr2_i_msa(spm, a1 = 1, a2 = 1))
#' }
#' @export
calculate_dr2_i_msa <- function(spm, a1, a2) {
  dr2_i <- dr2_msa(spm$dr2_ijm, spm$energy_data, a1, a2)
  tibble(i = seq_along(dr2_i), dr2_i = dr2_i)
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
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
#' @param a1 Stability selection strength (non-negative). `0` disables stability
#'   selection.
#' @param a2 Activity selection strength (non-negative). `0` disables activity
#'   selection.
#' @return A tibble with one row per site: the site index `i` and the predicted
#'   profile `lrmsd_i_msa`. Callers needing the PDB residue number join
#'   `spm$site_map` themselves (as [calculate_lrmsd_i_nested_models()] does).
#' @seealso [calculate_dr2_i_msa()] (the squared-divergence forward map it wraps);
#'   [calculate_lrmsd_i_nested_models()] (the four-variant analysis expansion).
#' @family model
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, seed = 1024)
#' head(calculate_lrmsd_i_msa(spm, a1 = 1, a2 = 1))
#' }
#' @export
calculate_lrmsd_i_msa <- function(spm, a1, a2) {
  lrmsd <- lrmsd_msa(spm$dr2_ijm, spm$energy_data, a1, a2)
  tibble(i = seq_along(lrmsd), lrmsd_i_msa = lrmsd)
}

#' Predicted mean-centred per-site log structural divergence at one selection strength
#'
#' Mean-centred counterpart of [calculate_lrmsd_i_msa()]: the per-site profile centred by
#' its own mean over the **full model support** (all model residues),
#' `nlrmsd_i = lrmsd_i - mean_S(lrmsd_i)`. This is the quantity the site ML fit is on -- the
#' likelihood centres both sides -- so it is the natural forward map for exploring a
#' scenario `(a1, a2)` on the same scale the fit reports, without fitting. It is also the
#' point profile that [predict_profiles()] bands: the two agree exactly (the
#' predictor's `_mean` column is this function's `nlrmsd_i_msa`).
#'
#' Centring is over all model residues, deliberately agnostic to which residues a dataset
#' happens to observe. The fit instead centres over the observation-matched residues, so a
#' fitted profile sits at a slightly different level by design; see [predict_profiles()].
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
#' @param a1 Stability selection strength (non-negative). `0` disables stability selection.
#' @param a2 Activity selection strength (non-negative). `0` disables activity selection.
#' @return A tibble with one row per site: the site index `i` and the mean-centred profile
#'   `nlrmsd_i_msa`.
#' @seealso [calculate_lrmsd_i_msa()] (the uncentred profile it centres);
#'   [predict_profiles()] (the same quantity with delta-method bands from a fit);
#'   [calculate_nlrmsd_n_msa()] (the mode-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, seed = 1024)
#' head(calculate_nlrmsd_i_msa(spm, a1 = 1, a2 = 1))
#' }
#' @export
calculate_nlrmsd_i_msa <- function(spm, a1, a2) {
  nlrmsd <- nlrmsd_msa(spm$dr2_ijm, spm$energy_data, a1, a2)
  tibble(i = seq_along(nlrmsd), nlrmsd_i_msa = nlrmsd)
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
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
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
#' spm <- generate_spm_data(znb_wt, seed = 1024)
#' head(calculate_dr2_n_msa(spm, a1 = 1, a2 = 1))
#' }
#' @export
calculate_dr2_n_msa <- function(spm, a1, a2) {
  dr2_n <- dr2_msa(spm$dr2_njm, spm$energy_data, a1, a2)
  tibble(n = seq_along(dr2_n), dr2_n = dr2_n)
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
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
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
#' spm <- generate_spm_data(znb_wt, seed = 1024)
#' head(calculate_lrmsd_n_msa(spm, a1 = 1, a2 = 1))
#' }
#' @export
calculate_lrmsd_n_msa <- function(spm, a1, a2) {
  lrmsd <- lrmsd_msa(spm$dr2_njm, spm$energy_data, a1, a2)
  tibble(n = seq_along(lrmsd), lrmsd_n_msa = lrmsd)
}

#' Predicted mean-centred per-mode log structural divergence at one selection strength
#'
#' Mode-indexed counterpart of [calculate_nlrmsd_i_msa()]: the per-mode profile centred by
#' its own mean over the full model support (all modes),
#' `nlrmsd_n = lrmsd_n - mean_S(lrmsd_n)`. This is the quantity the mode ML fit is on and the
#' point profile that [predict_profiles()] bands (the two agree exactly). Modes are
#' not residue-anchored, so there is no `pdb_site` column. See [calculate_nlrmsd_i_msa()] for
#' the support note (predict/scenario centres over all modes; the fit centres over the
#' observation-matched modes, so the two levels differ by design).
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
#' @param a1 Stability selection strength (non-negative). `0` disables stability selection.
#' @param a2 Activity selection strength (non-negative). `0` disables activity selection.
#' @return A tibble with one row per mode: the mode index `n` and the mean-centred profile
#'   `nlrmsd_n_msa`.
#' @seealso [calculate_lrmsd_n_msa()] (the uncentred profile it centres);
#'   [predict_profiles()] (the same quantity with delta-method bands from a fit);
#'   [calculate_nlrmsd_i_msa()] (the residue-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, seed = 1024)
#' head(calculate_nlrmsd_n_msa(spm, a1 = 1, a2 = 1))
#' }
#' @export
calculate_nlrmsd_n_msa <- function(spm, a1, a2) {
  nlrmsd <- nlrmsd_msa(spm$dr2_njm, spm$energy_data, a1, a2)
  tibble(n = seq_along(nlrmsd), nlrmsd_n_msa = nlrmsd)
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
#' Pass the four returned columns to `decompose_nested()` to split the
#' divergence into its mutation, stability, and activity contributions.
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
#' @param a1 Stability selection strength used for the MS and MSA variants
#'   (non-negative).
#' @param a2 Activity selection strength used for the MA and MSA variants
#'   (non-negative).
#' @return A tibble with one row per site: the site index `i`, the PDB residue
#'   number `pdb_site`, and the four divergence profiles `lrmsd_i_mm`,
#'   `lrmsd_i_ms`, `lrmsd_i_ma`, `lrmsd_i_msa`.
#' @seealso [calculate_dr2_i_msa()] (the single-variant prediction it calls);
#'   `decompose_nested()` (splits the four columns into contributions);
#'   [calculate_lrmsd_n_nested_models()] (the mode-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, seed = 1024)
#' head(calculate_lrmsd_i_nested_models(spm, a1 = 1, a2 = 1))
#' }
#' @export
calculate_lrmsd_i_nested_models <- function(spm, a1, a2) {
  v <- lrmsd_nested_models(spm$dr2_ijm, spm$energy_data, a1, a2)
  tibble(
    i           = seq_along(v$mm),
    lrmsd_i_mm  = v$mm,
    lrmsd_i_ms  = v$ms,
    lrmsd_i_ma  = v$ma,
    lrmsd_i_msa = v$msa
  ) %>%
    left_join(spm$site_map, by = "i") %>%
    select(i, pdb_site, everything())
}

#' Mean-centred per-site divergence profiles under all four model variants (MM, MS, MA, MSA)
#'
#' Mean-centred counterpart of [calculate_lrmsd_i_nested_models()]: the same four nested
#' variants, but each centred by **its own** mean over the full model support (all model
#' residues), `nlrmsd_i_v = lrmsd_i_v - mean_S(lrmsd_i_v)`. Each variant is centred
#' independently -- MM by the MM mean, MSA by the MSA mean -- matching how the fit and
#' [predict_decomposition()] treat them. This is the natural scenario forward map
#' for comparing the centred variants without fitting.
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
#' @param a1 Stability selection strength used for the MS and MSA variants (non-negative).
#' @param a2 Activity selection strength used for the MA and MSA variants (non-negative).
#' @return A tibble with one row per site: `i`, `pdb_site`, and the four mean-centred profiles
#'   `nlrmsd_i_mm`, `nlrmsd_i_ms`, `nlrmsd_i_ma`, `nlrmsd_i_msa`.
#' @seealso [calculate_lrmsd_i_nested_models()] (the uncentred variants it centres);
#'   [predict_decomposition()] (the same variants with delta-method bands);
#'   [calculate_nlrmsd_n_nested_models()] (the mode-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, seed = 1024)
#' head(calculate_nlrmsd_i_nested_models(spm, a1 = 1, a2 = 1))
#' }
#' @export
calculate_nlrmsd_i_nested_models <- function(spm, a1, a2) {
  v <- nlrmsd_nested_models(spm$dr2_ijm, spm$energy_data, a1, a2)
  tibble(
    i            = seq_along(v$mm),
    nlrmsd_i_mm  = v$mm,
    nlrmsd_i_ms  = v$ms,
    nlrmsd_i_ma  = v$ma,
    nlrmsd_i_msa = v$msa
  ) %>%
    left_join(spm$site_map, by = "i") %>%
    select(i, pdb_site, everything())
}

#' Per-site divergence decomposition at one (a1, a2)
#'
#' Forward decomposition on the site axis: at a single pair of selection strengths
#' `(a1, a2)` it splits the predicted structural-divergence profile into its
#' mutation, stability, and activity contributions (`phi_mut`, `phi_stab`,
#' `phi_act`). It packages the two-step recipe -- evaluate the four nested model
#' variants with [calculate_lrmsd_i_nested_models()], then apply the sequential
#' split `decompose_nested()` -- so a caller wanting phi at a point need
#' not run and thread the four variant columns by hand. The three contributions sum
#' exactly to the full-model profile `lrmsd_i_msa`.
#'
#' This is the single-point (forward) form. Propagating a fitted posterior over
#' `(a1, a2)` to phi with credible bands is done downstream by evaluating this
#' function across the posterior's nodes or draws.
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
#' @param a1 Stability selection strength (non-negative).
#' @param a2 Activity selection strength (non-negative).
#' @return A tibble with one row per site: the site index `i`, the PDB residue
#'   number `pdb_site`, and the three contributions `phi_mut`, `phi_stab`,
#'   `phi_act`.
#' @seealso [calculate_lrmsd_i_nested_models()] (the four variants it splits);
#'   `decompose_nested()` (the pure sequential split);
#'   [calculate_decomposition_n_msa()] (the mode-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, seed = 1024)
#' head(calculate_decomposition_i_msa(spm, a1 = 1, a2 = 1))
#' }
#' @export
calculate_decomposition_i_msa <- function(spm, a1, a2) {
  phi  <- lrmsd_msa_decomposition(spm$dr2_ijm, spm$energy_data, a1, a2)
  keys <- tibble(i = seq_along(phi$phi_mut)) %>%
    left_join(spm$site_map, by = "i") %>%
    select(i, pdb_site)
  dplyr::bind_cols(keys, tibble::as_tibble(phi))
}

#' Mean-centred per-site divergence decomposition at one (a1, a2)
#'
#' Mean-centred counterpart of [calculate_decomposition_i_msa()]: the three contributions of
#' the *centred* per-site profile `nlrmsd_i_msa`, each centred by its own mean over the full
#' model support and emitted as `nphi_mut`, `nphi_stab`, `nphi_act`. The decomposition is
#' defined on the centred quantity (the one the fit is on); the three columns sum exactly to
#' `nlrmsd_i_msa` (centring is linear). This is the point decomposition that
#' [predict_decomposition()] bands (the two agree exactly).
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
#' @param a1 Stability selection strength (non-negative).
#' @param a2 Activity selection strength (non-negative).
#' @return A tibble with one row per site: `i`, `pdb_site`, and the three mean-centred
#'   contributions `nphi_mut`, `nphi_stab`, `nphi_act`.
#' @seealso [calculate_decomposition_i_msa()] (the uncentred `phi_*` decomposition);
#'   [calculate_nlrmsd_i_msa()] (the centred profile these sum to);
#'   [predict_decomposition()] (the same contributions with delta-method
#'   bands); [calculate_nlrmsd_n_msa_decomposition()] (the mode-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, seed = 1024)
#' head(calculate_nlrmsd_i_msa_decomposition(spm, a1 = 1, a2 = 1))
#' }
#' @export
calculate_nlrmsd_i_msa_decomposition <- function(spm, a1, a2) {
  nphi <- nlrmsd_msa_decomposition(spm$dr2_ijm, spm$energy_data, a1, a2)
  tibble(i = seq_along(nphi$nphi_mut)) %>%
    left_join(spm$site_map, by = "i") %>%
    select(i, pdb_site) %>%
    dplyr::bind_cols(tibble::as_tibble(nphi))
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
#' `decompose_nested()` to split the divergence into its contributions.
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
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
#' spm <- generate_spm_data(znb_wt, seed = 1024)
#' head(calculate_lrmsd_n_nested_models(spm, a1 = 1, a2 = 1))
#' }
#' @export
calculate_lrmsd_n_nested_models <- function(spm, a1, a2) {
  v <- lrmsd_nested_models(spm$dr2_njm, spm$energy_data, a1, a2)
  tibble(
    n           = seq_along(v$mm),
    lrmsd_n_mm  = v$mm,
    lrmsd_n_ms  = v$ms,
    lrmsd_n_ma  = v$ma,
    lrmsd_n_msa = v$msa
  )
}

#' Mean-centred per-mode divergence profiles under all four model variants (MM, MS, MA, MSA)
#'
#' Mode-indexed counterpart of [calculate_nlrmsd_i_nested_models()]: the four nested variants
#' of the per-mode profile, each centred by its own mean over the full model support (all
#' modes). Modes are not residue-anchored, so there is no `pdb_site` column. Each variant is
#' centred independently, matching [predict_decomposition()].
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
#' @param a1 Stability selection strength used for the MS and MSA variants (non-negative).
#' @param a2 Activity selection strength used for the MA and MSA variants (non-negative).
#' @return A tibble with one row per mode: `n` and the four mean-centred profiles
#'   `nlrmsd_n_mm`, `nlrmsd_n_ms`, `nlrmsd_n_ma`, `nlrmsd_n_msa`.
#' @seealso [calculate_lrmsd_n_nested_models()] (the uncentred variants it centres);
#'   [predict_decomposition()] (the same variants with delta-method bands);
#'   [calculate_nlrmsd_i_nested_models()] (the residue-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, seed = 1024)
#' head(calculate_nlrmsd_n_nested_models(spm, a1 = 1, a2 = 1))
#' }
#' @export
calculate_nlrmsd_n_nested_models <- function(spm, a1, a2) {
  v <- nlrmsd_nested_models(spm$dr2_njm, spm$energy_data, a1, a2)
  tibble(
    n            = seq_along(v$mm),
    nlrmsd_n_mm  = v$mm,
    nlrmsd_n_ms  = v$ms,
    nlrmsd_n_ma  = v$ma,
    nlrmsd_n_msa = v$msa
  )
}

#' Per-mode divergence decomposition at one (a1, a2)
#'
#' Mode-indexed counterpart of [calculate_decomposition_i_msa()]: at a single pair
#' of selection strengths `(a1, a2)` it splits the predicted mode-divergence profile
#' into its mutation, stability, and activity contributions (`phi_mut`, `phi_stab`,
#' `phi_act`). It packages the two-step recipe -- evaluate the four nested model
#' variants with [calculate_lrmsd_n_nested_models()], then apply the sequential split
#' `decompose_nested()`. The three contributions sum exactly to the
#' full-model profile `lrmsd_n_msa`. Modes are not anchored to residues, so the
#' output has no `pdb_site` column.
#'
#' This is the single-point (forward) form. Propagating a fitted posterior over
#' `(a1, a2)` to phi with credible bands is done downstream by evaluating this
#' function across the posterior's nodes or draws.
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
#' @param a1 Stability selection strength (non-negative).
#' @param a2 Activity selection strength (non-negative).
#' @return A tibble with one row per mode: the mode index `n` and the three
#'   contributions `phi_mut`, `phi_stab`, `phi_act`.
#' @seealso [calculate_lrmsd_n_nested_models()] (the four variants it splits);
#'   `decompose_nested()` (the pure sequential split);
#'   [calculate_decomposition_i_msa()] (the residue-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, seed = 1024)
#' head(calculate_decomposition_n_msa(spm, a1 = 1, a2 = 1))
#' }
#' @export
calculate_decomposition_n_msa <- function(spm, a1, a2) {
  phi <- lrmsd_msa_decomposition(spm$dr2_njm, spm$energy_data, a1, a2)
  dplyr::bind_cols(tibble(n = seq_along(phi$phi_mut)), tibble::as_tibble(phi))
}

#' Mean-centred per-mode divergence decomposition at one (a1, a2)
#'
#' Mode-indexed counterpart of [calculate_nlrmsd_i_msa_decomposition()]: the three
#' contributions of the centred per-mode profile `nlrmsd_n_msa`, each centred by its own mean
#' over the full model support and emitted as `nphi_mut`, `nphi_stab`, `nphi_act`. The three
#' columns sum exactly to `nlrmsd_n_msa`. Modes are not residue-anchored, so there is no
#' `pdb_site` column. This is the point decomposition that
#' [predict_decomposition()] bands.
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
#' @param a1 Stability selection strength (non-negative).
#' @param a2 Activity selection strength (non-negative).
#' @return A tibble with one row per mode: `n` and the three mean-centred contributions
#'   `nphi_mut`, `nphi_stab`, `nphi_act`.
#' @seealso [calculate_decomposition_n_msa()] (the uncentred `phi_*` decomposition);
#'   [calculate_nlrmsd_n_msa()] (the centred profile these sum to);
#'   [predict_decomposition()] (the same contributions with delta-method
#'   bands); [calculate_nlrmsd_i_msa_decomposition()] (the residue-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, seed = 1024)
#' head(calculate_nlrmsd_n_msa_decomposition(spm, a1 = 1, a2 = 1))
#' }
#' @export
calculate_nlrmsd_n_msa_decomposition <- function(spm, a1, a2) {
  nphi <- nlrmsd_msa_decomposition(spm$dr2_njm, spm$energy_data, a1, a2)
  dplyr::bind_cols(tibble(n = seq_along(nphi$nphi_mut)), tibble::as_tibble(nphi))
}
