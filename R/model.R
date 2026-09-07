# MSA model forward map (divergence calculators at given (a1, a2))
# Functions to compute the model's predicted structural-divergence profiles

#' Fixation probability of a mutant under the MSA model
#'
#' The probability that a single-point mutant fixes, under stability selection of
#' strength `a1` and activity selection of strength `a2`.
#'
#' @details
#' The probability is
#' `p_fix = min(exp(-a1 * ddg), 1) * min(exp(-a2 * ddgact), 1)`, a property of one
#' mutant on its own: it depends only on that mutant's two energy changes and the
#' selection strengths, not on any ensemble of mutants.
#'
#' The function is vectorised. `ddg` and `ddgact` may be scalars, for one mutant, or
#' equal-length vectors, for many, and the result matches their shape.
#'
#' The value is not normalised. Turning fixation probabilities into averaging weights
#' over a particular ensemble of mutants is a separate step, which the profile
#' functions perform internally.
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
#' @seealso [calculate_profiles()] (the forward map built on these: it averages the
#'   per-mutant divergences, weighted by these fixation probabilities).
#' @examples
#' # One mutant, from its two energy changes:
#' pfix_msa(ddg = 1.2, ddgact = 0.4, a1 = 1, a2 = 1)
#'
#' # Vectorised over many mutants, and 0 switches a pressure off:
#' pfix_msa(ddg = c(0.5, 1.2, 3.0), ddgact = c(0.1, 0.4, 2.0), a1 = 1, a2 = 0)
#'
#' # A whole scan. `pdb_site_active` is required for the activity term: without it
#' # every `ddgact` is NA, and so is every fixation probability.
#' if (requireNamespace("bio3d", quietly = TRUE)) {
#'   ex  <- function(f) system.file("extdata", f, package = "msamodel")
#'   wt  <- penm::set_enm(bio3d::read.pdb(ex("1d6o_A.pdb")), node = "ca",
#'                               model = "ming_wall", d_max = 10.5, frustrated = FALSE)
#'   act <- read.csv(ex("1d6o_A_active_site.csv"))
#'   spm <- generate_spm(wt, pdb_site_active = act$pdb_site, ensemble = 1L)
#'   pfix_msa(spm$energy_data$ddg, spm$energy_data$ddgact, a1 = 1, a2 = 1)
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
#' `site_map` IS the key table: `(site, pdb_site)`, one row per site, already in
#' `dr2mat_site` column order, so it is bound on positionally rather than joined against a
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
#' `mode_map` is a single `mode` column, the index being the whole map, but it is a
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
#' itself stays user-facing: it says what happened, not what to go edit.
#'
#' @param metric The metric string that reached no branch.
#' @return A character message for `stop()`.
#' @noRd
unimplemented_metric_message <- function(metric) {
  paste0("metric \"", metric, "\" is not available for this function.")
}

# ---- the model layer public verbs: evaluate at a GIVEN (a1, a2), no fit ----------


#' Divergence profile at given selection strengths
#'
#' Calculates the structural divergence profile predicted by the MSA model at given
#' selection strengths, per residue and per normal mode.
#'
#' @details
#' The profile is `lrmsd_msa`, the divergence predicted by the full model, with both
#' selection pressures acting at the strengths given. `metric = "nlrmsd"` returns it
#' centred on its own mean instead, as `nlrmsd_msa`.
#'
#' The selection strengths are supplied here rather than estimated, so the result
#' carries no standard errors. For the profile at fitted strengths, with standard
#' errors, see [predict_profiles()].
#'
#' @param spm An `spm` object from [generate_spm()].
#' @param a1 Stability selection strength, non-negative. `0` switches stability
#'   selection off.
#' @param a2 Activity selection strength, non-negative. `0` switches activity
#'   selection off. It is on a different scale from `a1`, so the two numbers are not
#'   comparable to each other.
#' @param metric `"lrmsd"` for the profile as predicted, `"nlrmsd"` for the
#'   mean-centred one. Default `"lrmsd"`.
#' @return A list of two tibbles, `$site` (one row per residue, keyed by `site`, the
#'   column position in the scan, and `pdb_site`, the PDB residue number) and `$mode`
#'   (one row per normal mode, keyed by `mode`). Each carries one value column,
#'   `lrmsd_msa` or `nlrmsd_msa` according to `metric`.
#' @seealso [calculate_decomposition()] for the same profile split into the
#'   contributions of mutation, stability and activity selection;
#'   [predict_profiles()] for the profile at fitted selection strengths, with standard
#'   errors; [pfix_msa()] for the fixation probability `a1` and `a2` enter.
#' @family api
#' @examples
#' if (requireNamespace("bio3d", quietly = TRUE)) {
#'   ex  <- function(f) system.file("extdata", f, package = "msamodel")
#'   wt  <- penm::set_enm(bio3d::read.pdb(ex("1d6o_A.pdb")), node = "ca",
#'                        model = "ming_wall", d_max = 10.5, frustrated = FALSE)
#'   act <- read.csv(ex("1d6o_A_active_site.csv"))
#'   spm <- generate_spm(wt, pdb_site_active = act$pdb_site, ensemble = 1L)
#'
#'   calculate_profiles(spm, a1 = 2, a2 = 500)$site
#'
#'   # With both pressures off, the profile is mutation alone.
#'   calculate_profiles(spm, a1 = 0, a2 = 0)$site
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

#' Nested-model profiles and their mutation, stability and activity contributions
#'
#' Calculates the divergence profiles predicted by four nested models at given
#' selection strengths, and decomposes the full model's profile into mutation,
#' stability and activity contributions, per residue and per normal mode.
#'
#' @details
#' Setting a selection strength to zero switches that pressure off, giving four
#' models: MM `(0, 0)`, MS `(a1, 0)`, MA `(0, a2)` and MSA `(a1, a2)`. All four are
#' evaluated, as `lrmsd_mm`, `lrmsd_ms`, `lrmsd_ma` and `lrmsd_msa`. The full model's
#' profile is then split into increments along MM to MS to MSA:
#'
#' \preformatted{phi_mut  = lrmsd_mm
#' phi_stab = lrmsd_ms  - lrmsd_mm
#' phi_act  = lrmsd_msa - lrmsd_ms}
#'
#' These sum to `lrmsd_msa`. MA is outside the progression.
#'
#' `metric = "nlrmsd"` centres every column on its own mean, renaming the
#' contributions `nphi_mut`, `nphi_stab` and `nphi_act`.
#'
#' @param spm An `spm` object from [generate_spm()], carrying the per-mutant energies
#'   and squared displacements the profiles are averaged from.
#' @param a1 Stability selection strength, non-negative. It enters a mutant's fixation
#'   probability as `min(exp(-a1 * ddg), 1)`, so larger values suppress destabilising
#'   mutants more sharply, and `0` switches stability selection off.
#' @param a2 Activity selection strength, non-negative. It enters the same way, as
#'   `min(exp(-a2 * ddgact), 1)`, and `0` switches activity selection off. It is on a
#'   different scale from `a1`, so the two numbers are not comparable to each other.
#' @param metric `"lrmsd"` for the profiles as predicted, `"nlrmsd"` for the
#'   mean-centred ones. Default `"lrmsd"`.
#'
#' @return A list of two tibbles, `$site` (one row per residue, keyed by `site`, the
#'   column position in the scan, and `pdb_site`, the PDB residue number) and `$mode`
#'   (one row per normal mode, keyed by `mode`). Apart from those keys the two have the
#'   same columns.
#'
#'   With `metric = "lrmsd"`: `lrmsd_mm`, `lrmsd_ms`, `lrmsd_ma` and `lrmsd_msa`, the
#'   profile predicted by each of the four models, then `phi_mut`, `phi_stab` and
#'   `phi_act`, the three contributions, which sum to `lrmsd_msa`.
#'
#'   With `metric = "nlrmsd"`: the same seven columns, each mean-centred, named
#'   `nlrmsd_mm`, `nlrmsd_ms`, `nlrmsd_ma`, `nlrmsd_msa`, `nphi_mut`, `nphi_stab` and
#'   `nphi_act`. The three contributions sum to `nlrmsd_msa`.
#'
#' @seealso [calculate_profiles()] for the full model's profile alone, without the
#'   nested models and the contributions; [predict_decomposition()] for the same
#'   decomposition at fitted selection strengths, with standard errors;
#'   [fit_lrmsd_msa_site()] and [fit_lrmsd_msa_mode()] to estimate `a1` and `a2` from an
#'   observed profile. `vignette("msamodel-explore")` plots the nested models and the
#'   contributions in both representations.
#' @family api
#' @examples
#' if (requireNamespace("bio3d", quietly = TRUE)) {
#'   ex  <- function(f) system.file("extdata", f, package = "msamodel")
#'   wt  <- penm::set_enm(bio3d::read.pdb(ex("1d6o_A.pdb")), node = "ca",
#'                        model = "ming_wall", d_max = 10.5, frustrated = FALSE)
#'   act <- read.csv(ex("1d6o_A_active_site.csv"))
#'   spm <- generate_spm(wt, pdb_site_active = act$pdb_site, ensemble = 1L)
#'
#'   # Moderate stability selection, strong activity selection.
#'   d <- calculate_decomposition(spm, a1 = 2, a2 = 500, metric = "nlrmsd")
#'   d$site
#'
#'   # The residues where activity selection changes the profile most.
#'   head(d$site[order(d$site$nphi_act), c("pdb_site", "nphi_stab", "nphi_act")])
#'
#'   # The same decomposition in the per-mode representation.
#'   d$mode
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
