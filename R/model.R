#' MSA model forward map (divergence calculators at given (a1, a2))
#' Functions to compute the model's predicted structural-divergence profiles

#' Calculate dr2_i for a single selection parameter point using matrix operations
#'
#' @param spm_pp Output from [preprocess_spm()] (energy_data + the `dr2_ijm` matrix
#'   + site_map).
#' @param a1 Stability selection parameter value
#' @param a2 Activity selection parameter value
#' @return Tibble with columns i (site number) and dr2_i (weighted average)
#' @family model
#' @export
calculate_dr2_i_msa <- function(spm_pp, a1, a2) {
  energy_data <- spm_pp$energy_data
  dr2_ijm <- spm_pp$dr2_ijm

  # Fixation probabilities
  pstab_jm <- pmin(exp(-a1 * energy_data$ddg_jm), 1)
  pact_jm <- pmin(exp(-a2 * energy_data$ddgact_jm), 1)
  pfix_jm <- pstab_jm * pact_jm
  weights_jm <- pfix_jm / sum(pfix_jm)

  # Weighted average over the mutant (j,m) axis: (dr2_ijm) -> (dr2_i)
  dr2_i <- colSums(dr2_ijm * weights_jm)

  tibble(i = as.integer(colnames(dr2_ijm)), dr2_i = dr2_i)
}

#' Calculate dr2_n (mode-form structural divergence) for a selection-parameter point
#'
#' Mode counterpart of [calculate_dr2_i_msa()]. The reweighting is identical --
#' the fixation weights live only on the mutant axis, so the same
#' `colSums(value_matrix * weights_jm)` collapses the mutant rows whether the
#' columns are sites or modes. Predict-only: there is no observed mode profile to
#' fit, so this has no log-likelihood counterpart.
#'
#' @param spm_pp Output from [preprocess_spm_mode()] (energy_data + the `dr2_njm`
#'   matrix).
#' @param a1 Stability selection parameter value
#' @param a2 Activity selection parameter value
#' @return Tibble with columns `n` (mode index) and `dr2_n` (weighted average
#'   per-mode squared divergence)
#' @family model
#' @export
calculate_dr2_n_msa <- function(spm_pp, a1, a2) {
  energy_data <- spm_pp$energy_data
  dr2_njm <- spm_pp$dr2_njm

  # Fixation probabilities (same axis-agnostic weights as the site form)
  pstab_jm <- pmin(exp(-a1 * energy_data$ddg_jm), 1)
  pact_jm <- pmin(exp(-a2 * energy_data$ddgact_jm), 1)
  pfix_jm <- pstab_jm * pact_jm
  weights_jm <- pfix_jm / sum(pfix_jm)

  # Weighted average over the mutant (j,m) axis: (dr2_njm) -> (dr2_n)
  dr2_n <- colSums(dr2_njm * weights_jm)

  tibble(n = as.integer(colnames(dr2_njm)), dr2_n = dr2_n)
}
