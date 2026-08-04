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

