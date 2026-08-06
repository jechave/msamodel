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
#' @param ddg Stability free-energy change(s) of the mutant(s) (`ddg_jm`). Scalar or
#'   vector.
#' @param ddgact Activity free-energy change(s) of the mutant(s) (`ddgact_jm`), the
#'   same length as `ddg`.
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
#' @return A numeric vector of averaging weights, one per mutant, summing to one.
#' @family model
#' @noRd
weights_jm <- function(energy_data, a1, a2) {
  pfix_jm <- pfix_msa(energy_data$ddg_jm, energy_data$ddgact_jm, a1, a2)
  pfix_jm / sum(pfix_jm)
}

#' Axis-blind per-response structural divergence at one (a1, a2)
#'
#' Weights each mutant by its MSA fixation probability and averages the per-response
#' squared displacements over mutants. The forward-map core the leaf verbs are built on.
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
#' The four variants MM `(0,0)`, MS `(a1,0)`, MA `(0,a2)`, MSA `(a1,a2)`, each an
#' `lrmsd` vector.
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
#' Each nested variant centred by its own mean.
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
#' Evaluate the four nested variants, then apply the sequential split
#' `decompose_nested()`. Named `lrmsd_msa_decomposition` so it pairs cleanly with
#' `nlrmsd_msa_decomposition`.
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
#' Each phi contribution centred by its own mean.
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

# ---- assembly helpers: bare vectors -> the keyed tibble the verbs return ----------

#' Attach the axis key (and, for site, pdb_site) to bare value vectors
#'
#' The forward-map primitives return bare per-response vectors; this wraps a NAMED list
#' of equal-length bare vectors into the keyed tibble the verbs return. The list names
#' become the value-column names verbatim -- the verb owns naming (including any `_se`
#' siblings); `key_profile` only prepends the axis key and, for the site axis, joins
#' `pdb_site` via `spm$site_map`.
#'
#' On the site axis `site_map` IS the key table -- `(site, pdb_site)`, one row per site,
#' in `dr2_ijm` column order -- so it is bound on positionally rather than joined against
#' a manufactured index. The row-count equality that makes that valid is asserted, not
#' assumed: a positional bind fails loud on a mismatch where a join would have silently
#' filled `pdb_site` with `NA`.
#'
#' @param cols A NAMED list of equal-length bare numeric vectors. Names become columns.
#' @param spm The `spm` object (for `site_map` on the site axis).
#' @param axis `"site"` (keys `site`, `pdb_site`) or `"mode"` (key `mode`).
#' @return A tibble: site -> `site, pdb_site, <cols...>`; mode -> `mode, <cols...>`.
#' @noRd
key_profile <- function(cols, spm, axis) {
  n <- length(cols[[1]])
  body <- tibble::as_tibble(cols)
  if (axis == "site") {
    if (nrow(spm$site_map) != n) {
      stop("site_map has ", nrow(spm$site_map), " rows but the profile has ", n,
           " sites; the spm object is inconsistent.")
    }
    dplyr::bind_cols(spm$site_map, body)
  } else {
    dplyr::bind_cols(tibble(mode = seq_len(n)), body)
  }
}

#' Per-axis dispatch table: the divergence matrix for one branch
#'
#' The single place the site/mode branch is chosen: `dr2_ijm` + `"site"` or
#' `dr2_njm` + `"mode"`. The verbs iterate over the two branches and call the axis-blind
#' primitives with the branch's matrix.
#'
#' Carries no name tag: both branches emit the SAME value-column vocabulary
#' (`lrmsd_msa`, `lrmsd_mm`, ...), so only the key column differs, and `key_profile()`
#' derives that from `axis`.
#'
#' @param spm The `spm` object.
#' @return A named list of two branches, each `list(dr2_mat =, axis =)`.
#' @noRd
axis_branches <- function(spm) {
  list(site = list(dr2_mat = spm$dr2_ijm, axis = "site"),
       mode = list(dr2_mat = spm$dr2_njm, axis = "mode"))
}

# ---- the model layer public verbs: evaluate at a GIVEN (a1, a2), no fit ----------


#' Predicted divergence profiles at one selection strength (both axes)
#'
#' The model's forward per-response log structural-divergence profile at a single pair
#' of selection strengths `(a1, a2)`, on **both** response axes at once. `which` selects
#' the quantity: `"lrmsd"` (the absolute profile `log(sqrt(dr2))`) or `"nlrmsd"` (the
#' mean-centred profile the fit is on). Returns point values only -- no error bands
#' (there is no fit and hence no parameter covariance); for bands from a fit use
#' [predict_profiles()].
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
#' @param a1 Stability selection strength (non-negative). `0` disables it.
#' @param a2 Activity selection strength (non-negative). `0` disables it.
#' @param which `"lrmsd"` (absolute) or `"nlrmsd"` (mean-centred). Default `"lrmsd"`.
#' @return A list with two tibbles. `$site`: `site`, `pdb_site`, and the profile column
#'   (`lrmsd_msa` or `nlrmsd_msa`). `$mode`: `mode` and the same profile column. Both
#'   branches use identical value-column names; only the key column differs.
#' @seealso [predict_profiles()] (the same profiles with error bands, from a fit);
#'   [calculate_decomposition()] (the profile split into contributions).
#' @family api
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' calculate_profiles(spm, a1 = 1, a2 = 1, which = "nlrmsd")$site
#' }
#' @export
calculate_profiles <- function(spm, a1, a2, which = c("lrmsd", "nlrmsd")) {
  which <- match.arg(which)
  fwd   <- if (which == "nlrmsd") nlrmsd_msa else lrmsd_msa
  out <- lapply(axis_branches(spm), function(b) {
    v    <- fwd(b$dr2_mat, spm$energy_data, a1, a2)
    name <- paste0(which, "_msa")
    key_profile(stats::setNames(list(v), name), spm, b$axis)
  })
  out
}

# ---- calculate_decomposition -----------------------------------------------------

#' Divergence decomposition at one selection strength (both axes)
#'
#' The nested-model profiles AND the three sequential contributions of the divergence
#' profile at a single `(a1, a2)`, on **both** response axes. `which` selects the family
#' in lockstep: `"lrmsd"` gives the uncentred nested models (`lrmsd_mm`...) and the
#' uncentred contributions (`phi_mut`, `phi_stab`, `phi_act`); `"nlrmsd"` gives the
#' mean-centred nested models (`nlrmsd_mm`...) and the centred contributions
#' (`nphi_mut`, `nphi_stab`, `nphi_act`). Point values only; for bands from a fit use
#' [predict_decomposition()].
#'
#' The three contributions sum exactly to the full-model (`msa`) profile.
#'
#' @param spm A single-point-mutation `spm` object from [generate_spm_data()].
#' @param a1 Stability selection strength (non-negative).
#' @param a2 Activity selection strength (non-negative).
#' @param which `"lrmsd"` (absolute) or `"nlrmsd"` (mean-centred). Default `"lrmsd"`.
#' @return A list with two tibbles (`$site`, `$mode`). Each holds the axis key
#'   (`site`, `pdb_site` for site; `mode` for mode), the four nested-model columns, and
#'   the three contribution columns, all in the family selected by `which`. Both branches
#'   use identical value-column names.
#' @seealso [predict_decomposition()] (the same, with error bands from a fit);
#'   [calculate_profiles()] (the profile these contributions sum to).
#' @family api
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, pdb_site_active = c(99,101,103,162,181,184,193,223), seed = 1024)
#' calculate_decomposition(spm, a1 = 1, a2 = 1, which = "nlrmsd")$site
#' }
#' @export
calculate_decomposition <- function(spm, a1, a2, which = c("lrmsd", "nlrmsd")) {
  which <- match.arg(which)
  centred     <- which == "nlrmsd"
  nested_fwd  <- if (centred) nlrmsd_nested_models else lrmsd_nested_models
  decomp_fwd  <- if (centred) nlrmsd_msa_decomposition else lrmsd_msa_decomposition
  phi_prefix  <- if (centred) "nphi" else "phi"

  out <- lapply(axis_branches(spm), function(b) {
    nm  <- nested_fwd(b$dr2_mat, spm$energy_data, a1, a2)          # list mm/ms/ma/msa
    phi <- decomp_fwd(b$dr2_mat, spm$energy_data, a1, a2)          # list *_mut/stab/act
    nested_cols <- stats::setNames(
      nm[c("mm", "ms", "ma", "msa")],
      paste0(which, "_", c("mm", "ms", "ma", "msa")))
    # decomp_fwd already names phi_mut/... (lrmsd) or nphi_mut/... (nlrmsd) -- take as is.
    key_profile(c(nested_cols, phi), spm, b$axis)
  })
  out
}
