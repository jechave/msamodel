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
#' primitive the ensemble averaging weights (`weights_jm()`) are built from.
#'
#' Pure and vectorised: `ddg` and `ddgact` may be scalars (one mutant) or
#' equal-length vectors (many mutants), and the result matches their shape. This is
#' *not* normalised -- turning fixation probabilities into averaging weights over a
#' particular ensemble of mutants is a separate, ensemble-specific step.
#'
#' @param ddg Stability free-energy change(s) of the mutant(s), as carried in an
#'   `spm` object's `energy_data$ddg`. Scalar or vector.
#' @param ddgact Activity free-energy change(s) of the mutant(s)
#'   (`energy_data$ddgact`), the same length as `ddg`.
#' @param a1 Stability selection strength (non-negative). `0` disables stability
#'   selection.
#' @param a2 Activity selection strength (non-negative). `0` disables activity
#'   selection.
#' @return A numeric vector of fixation probabilities, the same length as `ddg`.
#' @seealso [calculate_profiles()] (the forward map built on these; it averages the
#'   per-mutant divergences with the normalised `weights_jm()`).
#' @family model
#' @examples
#' \dontrun{
#' # One mutant:
#' pfix_msa(ddg = 1.2, ddgact = 0.4, a1 = 1, a2 = 1)
#'
#' # A whole scan:
#' ex  <- function(f) system.file("extdata", f, package = "msamodel")
#' wt  <- set_enm(bio3d::read.pdb(ex("1znb_A.pdb")), node = "ca",
#'                       model = "ming_wall", d_max = 10.5, frustrated = FALSE)
#' spm <- generate_spm(wt, ensemble = 1L)
#' pfix_msa(spm$energy_data$ddg, spm$energy_data$ddgact, a1 = 1, a2 = 1)
#' }
#' @export
pfix_msa <- function(ddg, ddgact, a1, a2) {
  pstab <- pmin(exp(-a1 * ddg), 1)
  pact  <- pmin(exp(-a2 * ddgact), 1)
  pstab * pact
}

# ---- axis-blind forward-map primitives (internal) ---------------------------------
# The site (_i_/dr2mat_site) and mode (_n_/dr2mat_mode) forward maps are the same math on a
# different response-mutation matrix. These primitives hold that shared math: they take
# a bare [mutant x response] matrix `dr2mat` (either dr2mat_site or dr2mat_mode) plus the
# `energy_data` tibble, and return bare vectors (or a named list of bare vectors) in
# matrix-column order -- no index, no tibble, no site_map. Alignment is by column
# POSITION (the dr2 matrices carry no column names by construction); the exported verbs
# below attach the axis key at the boundary, via prepend_site_key / prepend_mode_key.

#' SPM-ensemble averaging weights from an energy table
#'
#' Turns the per-mutant MSA fixation probabilities into the normalised averaging weights
#' `weights_jm = pfix_jm / sum(pfix_jm)`, one per mutant `(j, m)`, summing to one. Reads
#' the per-mutant energies directly from an `energy_data` tibble, so the forward-map
#' primitives obtain their weights without the `spm` object.
#'
#' @param energy_data A tibble carrying per-mutant `ddg` and `ddgact` columns.
#' @param a1 Stability selection strength (non-negative).
#' @param a2 Activity selection strength (non-negative).
#' @return A numeric vector of averaging weights, one per mutant, summing to one.
#' @family model
#' @noRd
weights_jm <- function(energy_data, a1, a2) {
  pfix_jm <- pfix_msa(energy_data$ddg, energy_data$ddgact, a1, a2)
  pfix_jm / sum(pfix_jm)
}

#' Axis-blind per-response structural divergence at one (a1, a2)
#'
#' Weights each mutant by its MSA fixation probability and averages the per-response
#' squared displacements over mutants. The forward-map core the leaf verbs are built on.
#'
#' @param dr2mat A `[mutant x response]` divergence matrix (`dr2mat_site` or `dr2mat_mode`).
#' @param energy_data The per-mutant energy tibble (for the weights).
#' @param a1,a2 Selection strengths.
#' @return A numeric vector of per-response `dr2`, in column order.
#' @family model
#' @noRd
dr2_msa <- function(dr2mat, energy_data, a1, a2) {
  w <- weights_jm(energy_data, a1, a2)
  colSums(dr2mat * w)
}

#' Axis-blind per-response log structural divergence at one (a1, a2)
#'
#' `lrmsd = log(sqrt(dr2))`. Sole owner of the `dr2 -> lrmsd` transform.
#'
#' @inheritParams dr2_msa
#' @return A numeric vector of per-response `lrmsd`, in column order.
#' @family model
#' @noRd
lrmsd_msa <- function(dr2mat, energy_data, a1, a2) {
  log(sqrt(dr2_msa(dr2mat, energy_data, a1, a2)))
}

#' Axis-blind mean-centred per-response log divergence at one (a1, a2)
#'
#' `nlrmsd = lrmsd - mean(lrmsd)` over the full response support.
#'
#' @inheritParams dr2_msa
#' @return A numeric vector of per-response `nlrmsd`, in column order.
#' @family model
#' @noRd
nlrmsd_msa <- function(dr2mat, energy_data, a1, a2) {
  lrmsd <- lrmsd_msa(dr2mat, energy_data, a1, a2)
  lrmsd - mean(lrmsd)
}

#' Axis-blind four nested-model lrmsd profiles at one (a1, a2)
#'
#' The four variants MM `(0,0)`, MS `(a1,0)`, MA `(0,a2)`, MSA `(a1,a2)`, each an
#' `lrmsd` vector.
#'
#' @inheritParams dr2_msa
#' @return A named list `mm`, `ms`, `ma`, `msa`, each an `lrmsd` vector in column order.
#' @family model
#' @noRd
lrmsd_nested_models <- function(dr2mat, energy_data, a1, a2) {
  lrmsd <- function(p1, p2) lrmsd_msa(dr2mat, energy_data, p1, p2)
  list(mm  = lrmsd(0,  0),
       ms  = lrmsd(a1, 0),
       ma  = lrmsd(0,  a2),
       msa = lrmsd(a1, a2))
}

#' Axis-blind four mean-centred nested-model nlrmsd profiles at one (a1, a2)
#'
#' Each nested variant centred by its own mean.
#'
#' @inheritParams dr2_msa
#' @return A named list `mm`, `ms`, `ma`, `msa`, each a centred vector in column order.
#' @family model
#' @noRd
nlrmsd_nested_models <- function(dr2mat, energy_data, a1, a2) {
  v <- lrmsd_nested_models(dr2mat, energy_data, a1, a2)
  lapply(v, function(x) x - mean(x))
}

#' Axis-blind divergence decomposition at one (a1, a2)
#'
#' Evaluate the four nested variants, then apply the sequential split
#' `decompose_nested()`. Named `lrmsd_msa_decomposition` so it pairs cleanly with
#' `nlrmsd_msa_decomposition`.
#'
#' @inheritParams dr2_msa
#' @return A list `phi_mut`, `phi_stab`, `phi_act`, each a vector in column order.
#' @family model
#' @noRd
lrmsd_msa_decomposition <- function(dr2mat, energy_data, a1, a2) {
  v <- lrmsd_nested_models(dr2mat, energy_data, a1, a2)
  decompose_nested(v$mm, v$ms, v$ma, v$msa)
}

#' Axis-blind mean-centred divergence decomposition at one (a1, a2)
#'
#' Each phi contribution centred by its own mean.
#'
#' @inheritParams dr2_msa
#' @return A named list `nphi_mut`, `nphi_stab`, `nphi_act`, each a centred vector in column order.
#' @family model
#' @noRd
nlrmsd_msa_decomposition <- function(dr2mat, energy_data, a1, a2) {
  phi <- lrmsd_msa_decomposition(dr2mat, energy_data, a1, a2)
  list(nphi_mut  = phi$phi_mut  - mean(phi$phi_mut),
       nphi_stab = phi$phi_stab - mean(phi$phi_stab),
       nphi_act  = phi$phi_act  - mean(phi$phi_act))
}

# ---- assembly helpers: bare vectors -> the keyed tibble the verbs return ----------

#' Prepend the site key to a site-axis value tibble
#'
#' `site_map` IS the key table -- `(site, pdb_site)`, one row per site, already in
#' `dr2mat_site` column order -- so it is bound on positionally rather than joined against a
#' manufactured index. The row-count equality that makes that valid is asserted, not
#' assumed: a positional bind fails loud on a mismatch where a join would have silently
#' filled `pdb_site` with `NA`.
#'
#' @param site_map The `(site, pdb_site)` key tibble, i.e. `spm$site_map`.
#' @param body A tibble of value columns, one row per site, in `dr2mat_site` column order.
#' @return `body` with `site` and `pdb_site` prepended.
#' @noRd
prepend_site_key <- function(site_map, body) {
  if (nrow(site_map) != nrow(body)) {
    stop("site_map has ", nrow(site_map), " rows but the profile has ", nrow(body),
         " sites; the spm object is inconsistent.")
  }
  dplyr::bind_cols(site_map, body)
}

#' Prepend the mode key to a mode-axis value tibble
#'
#' The mode-axis counterpart of `prepend_site_key()`. Modes are not residue-anchored, so
#' `mode_map` is a single `mode` column -- the index is the whole map -- but it is a
#' STORED key bound on positionally, and its row count is asserted for the same reason.
#'
#' @param mode_map The `(mode)` key tibble, i.e. `spm$mode_map`.
#' @param body A tibble of value columns, one row per mode, in `dr2mat_mode` column order.
#' @return `body` with `mode` prepended.
#' @noRd
prepend_mode_key <- function(mode_map, body) {
  if (nrow(mode_map) != nrow(body)) {
    stop("mode_map has ", nrow(mode_map), " rows but the profile has ", nrow(body),
         " modes; the spm object is inconsistent.")
  }
  dplyr::bind_cols(mode_map, body)
}

#' Error message for a metric that `match.arg` accepted but no branch implements
#'
#' Error message for a metric that `match.arg` accepted but no branch implements
#'
#' The metric-dispatching verbs branch explicitly on each metric they implement and end
#' in a `stop()` rather than an `else` that means "lrmsd". `match.arg` already rejects
#' misspellings, so this cannot fire through the public API today. It exists for the
#' case `match.arg` cannot catch: a metric ADDED to a verb's formals (`drmsf`, `dnh`,
#' ...) whose branch was never written. Without it the unwritten metric would silently
#' return lrmsd numbers under its own name.
#'
#' TO THE MAINTAINER, if you are reading this from a traceback: add the missing branch
#' to the verb that raised it, or drop the metric from that verb's formals. The message
#' itself stays user-facing -- it says what happened, not what to go edit.
#'
#' @param metric The metric string that reached no branch.
#' @return A character message for `stop()`.
#' @noRd
unimplemented_metric_message <- function(metric) {
  paste0("metric \"", metric, "\" is not available for this function.")
}

# ---- the model layer public verbs: evaluate at a GIVEN (a1, a2), no fit ----------


#' Predicted divergence profiles at one selection strength (both axes)
#'
#' The model's forward per-response log structural-divergence profile at a single pair
#' of selection strengths `(a1, a2)`, on **both** response axes at once. `metric` selects
#' the quantity: `"lrmsd"` (the absolute profile `log(sqrt(dr2))`) or `"nlrmsd"` (the
#' mean-centred profile the fit is on). Returns point values only -- no error bands
#' (there is no fit and hence no parameter covariance); for bands from a fit use
#' [predict_profiles()].
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm()].
#' @param a1 Stability selection strength (non-negative). `0` disables it.
#' @param a2 Activity selection strength (non-negative). `0` disables it.
#' @param metric `"lrmsd"` (absolute) or `"nlrmsd"` (mean-centred). Default `"lrmsd"`.
#' @return A list with two tibbles. `$site`: `site`, `pdb_site`, and the profile column
#'   (`lrmsd_msa` or `nlrmsd_msa`). `$mode`: `mode` and the same profile column. Both
#'   branches use identical value-column names; only the key column differs.
#' @seealso [predict_profiles()] (the same profiles with error bands, from a fit);
#'   [calculate_decomposition()] (the profile split into contributions).
#' @family api
#' @examples
#' \dontrun{
#' ex  <- function(f) system.file("extdata", f, package = "msamodel")
#' wt  <- set_enm(bio3d::read.pdb(ex("1znb_A.pdb")), node = "ca",
#'                       model = "ming_wall", d_max = 10.5, frustrated = FALSE)
#' act <- readr::read_csv(ex("znb_active_site.csv"))
#' spm <- generate_spm(wt, pdb_site_active = act$pdb_site, ensemble = 1L)
#' calculate_profiles(spm, a1 = 1, a2 = 1, metric = "nlrmsd")$site
#' }
#' @export
calculate_profiles <- function(spm, a1, a2, metric = c("lrmsd", "nlrmsd")) {
  metric <- match.arg(metric)

  # --- site axis (responses are residues; dr2mat_site)
  if (metric == "nlrmsd") {
    site_profile <- tibble(nlrmsd_msa = nlrmsd_msa(spm$dr2mat_site, spm$energy_data, a1, a2))
  } else if (metric == "lrmsd") {
    site_profile <- tibble(lrmsd_msa = lrmsd_msa(spm$dr2mat_site, spm$energy_data, a1, a2))
  } else {
    stop(unimplemented_metric_message(metric))
  }
  site_profile <- prepend_site_key(spm$site_map, site_profile)

  # --- mode axis (responses are normal modes; dr2mat_mode)
  if (metric == "nlrmsd") {
    mode_profile <- tibble(nlrmsd_msa = nlrmsd_msa(spm$dr2mat_mode, spm$energy_data, a1, a2))
  } else if (metric == "lrmsd") {
    mode_profile <- tibble(lrmsd_msa = lrmsd_msa(spm$dr2mat_mode, spm$energy_data, a1, a2))
  } else {
    stop(unimplemented_metric_message(metric))
  }
  mode_profile <- prepend_mode_key(spm$mode_map, mode_profile)

  list(site = site_profile, mode = mode_profile)
}

# ---- calculate_decomposition -----------------------------------------------------

#' Divergence decomposition at one selection strength (both axes)
#'
#' The nested-model profiles AND the three sequential contributions of the divergence
#' profile at a single `(a1, a2)`, on **both** response axes. `metric` applies to every
#' returned column at once -- you never get some columns absolute and others centred.
#' `"lrmsd"` returns the absolute nested models (`lrmsd_mm`...) with the absolute
#' contributions (`phi_mut`, `phi_stab`, `phi_act`); `"nlrmsd"` returns the mean-centred
#' nested models (`nlrmsd_mm`...) with the centred contributions (`nphi_mut`,
#' `nphi_stab`, `nphi_act`). Point values only; for bands from a fit use
#' [predict_decomposition()].
#'
#' The three contributions sum exactly to the full-model (`msa`) profile.
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm()].
#' @param a1 Stability selection strength (non-negative).
#' @param a2 Activity selection strength (non-negative).
#' @param metric `"lrmsd"` (absolute) or `"nlrmsd"` (mean-centred). Default `"lrmsd"`.
#' @return A list with two tibbles (`$site`, `$mode`). Each holds the axis key
#'   (`site`, `pdb_site` for site; `mode` for mode), the four nested-model columns, and
#'   the three contribution columns, all on the `metric` you asked for. Both branches
#'   use identical value-column names.
#' @seealso [predict_decomposition()] (the same, with error bands from a fit);
#'   [calculate_profiles()] (the profile these contributions sum to).
#' @family api
#' @examples
#' \dontrun{
#' ex  <- function(f) system.file("extdata", f, package = "msamodel")
#' wt  <- set_enm(bio3d::read.pdb(ex("1znb_A.pdb")), node = "ca",
#'                       model = "ming_wall", d_max = 10.5, frustrated = FALSE)
#' act <- readr::read_csv(ex("znb_active_site.csv"))
#' spm <- generate_spm(wt, pdb_site_active = act$pdb_site, ensemble = 1L)
#' calculate_decomposition(spm, a1 = 1, a2 = 1, metric = "nlrmsd")$site
#' }
#' @export
calculate_decomposition <- function(spm, a1, a2, metric = c("lrmsd", "nlrmsd")) {
  metric <- match.arg(metric)

  # --- site axis (responses are residues; dr2mat_site)
  if (metric == "nlrmsd") {
    nested        <- nlrmsd_nested_models(spm$dr2mat_site, spm$energy_data, a1, a2)
    decomposition <- nlrmsd_msa_decomposition(spm$dr2mat_site, spm$energy_data, a1, a2)
    site_decomposition <- tibble(
      nlrmsd_mm  = nested$mm,
      nlrmsd_ms  = nested$ms,
      nlrmsd_ma  = nested$ma,
      nlrmsd_msa = nested$msa,
      nphi_mut   = decomposition$nphi_mut,
      nphi_stab  = decomposition$nphi_stab,
      nphi_act   = decomposition$nphi_act)
  } else if (metric == "lrmsd") {
    nested        <- lrmsd_nested_models(spm$dr2mat_site, spm$energy_data, a1, a2)
    decomposition <- lrmsd_msa_decomposition(spm$dr2mat_site, spm$energy_data, a1, a2)
    site_decomposition <- tibble(
      lrmsd_mm  = nested$mm,
      lrmsd_ms  = nested$ms,
      lrmsd_ma  = nested$ma,
      lrmsd_msa = nested$msa,
      phi_mut   = decomposition$phi_mut,
      phi_stab  = decomposition$phi_stab,
      phi_act   = decomposition$phi_act)
  } else {
    stop(unimplemented_metric_message(metric))
  }
  site_decomposition <- prepend_site_key(spm$site_map, site_decomposition)

  # --- mode axis (responses are normal modes; dr2mat_mode)
  if (metric == "nlrmsd") {
    nested        <- nlrmsd_nested_models(spm$dr2mat_mode, spm$energy_data, a1, a2)
    decomposition <- nlrmsd_msa_decomposition(spm$dr2mat_mode, spm$energy_data, a1, a2)
    mode_decomposition <- tibble(
      nlrmsd_mm  = nested$mm,
      nlrmsd_ms  = nested$ms,
      nlrmsd_ma  = nested$ma,
      nlrmsd_msa = nested$msa,
      nphi_mut   = decomposition$nphi_mut,
      nphi_stab  = decomposition$nphi_stab,
      nphi_act   = decomposition$nphi_act)
  } else if (metric == "lrmsd") {
    nested        <- lrmsd_nested_models(spm$dr2mat_mode, spm$energy_data, a1, a2)
    decomposition <- lrmsd_msa_decomposition(spm$dr2mat_mode, spm$energy_data, a1, a2)
    mode_decomposition <- tibble(
      lrmsd_mm  = nested$mm,
      lrmsd_ms  = nested$ms,
      lrmsd_ma  = nested$ma,
      lrmsd_msa = nested$msa,
      phi_mut   = decomposition$phi_mut,
      phi_stab  = decomposition$phi_stab,
      phi_act   = decomposition$phi_act)
  } else {
    stop(unimplemented_metric_message(metric))
  }
  mode_decomposition <- prepend_mode_key(spm$mode_map, mode_decomposition)

  list(site = site_decomposition, mode = mode_decomposition)
}
