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

#' Per-site lrmsd profiles for the four nested MSA model variants
#'
#' Evaluates the forward map [calculate_dr2_i_msa()] at the four nested
#' selection-parameter points and returns their per-site lrmsd profiles, where
#' `lrmsd = log(sqrt(dr2_i))`:
#'
#' | variant | (a1, a2) | meaning |
#' |---------|----------|---------|
#' | MM  | (0, 0)   | mutation only -- no selection |
#' | MS  | (a1, 0)  | mutation + stability selection |
#' | MA  | (0, a2)  | mutation + activity selection |
#' | MSA | (a1, a2) | full model -- both selections |
#'
#' This is the single source of truth for the four-variant recipe: both the
#' MCMC path ([calculate_prediction_samples()], per posterior sample) and the
#' fixed-(a1, a2) analysis use it. Feed the four columns to
#' [calculate_msa_decomposition()] for the phi decomposition.
#'
#' @param spm_pp Output from [preprocess_spm()] (energy_data + `dr2_ijm` + site_map).
#' @param a1 Stability selection parameter value.
#' @param a2 Activity selection parameter value.
#' @return Tibble with `i`, `pdb_site`, and the four nested-model lrmsd columns
#'   `lrmsd_i_mm`, `lrmsd_i_ms`, `lrmsd_i_ma`, `lrmsd_i_msa`.
#' @family model
#' @export
calculate_lrmsd_i_nested_models <- function(spm_pp, a1, a2) {
  lrmsd <- function(p1, p2) log(sqrt(calculate_dr2_i_msa(spm_pp, p1, p2)$dr2_i))

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

#' Per-mode lrmsd profiles for the four nested MSA model variants
#'
#' Mode counterpart of [calculate_lrmsd_i_nested_models()]. Evaluates the mode
#' forward map [calculate_dr2_n_msa()] at the four nested selection-parameter
#' points and returns their per-mode lrmsd profiles, where
#' `lrmsd = log(sqrt(dr2_n))`:
#'
#' | variant | (a1, a2) | meaning |
#' |---------|----------|---------|
#' | MM  | (0, 0)   | mutation only -- no selection |
#' | MS  | (a1, 0)  | mutation + stability selection |
#' | MA  | (0, a2)  | mutation + activity selection |
#' | MSA | (a1, a2) | full model -- both selections |
#'
#' Single source of truth for the mode four-variant recipe. Feed the four columns
#' to [calculate_msa_decomposition()] for the phi decomposition. Modes are not
#' anchored to residues, so there is no `pdb_site` (unlike the site form).
#'
#' @param spm_pp Output from [preprocess_spm_mode()] (energy_data + `dr2_njm`).
#' @param a1 Stability selection parameter value.
#' @param a2 Activity selection parameter value.
#' @return Tibble with `n` (mode index) and the four nested-model lrmsd columns
#'   `lrmsd_n_mm`, `lrmsd_n_ms`, `lrmsd_n_ma`, `lrmsd_n_msa`.
#' @family model
#' @export
calculate_lrmsd_n_nested_models <- function(spm_pp, a1, a2) {
  lrmsd <- function(p1, p2) log(sqrt(calculate_dr2_n_msa(spm_pp, p1, p2)$dr2_n))

  tibble(
    n           = as.integer(colnames(spm_pp$dr2_njm)),
    lrmsd_n_mm  = lrmsd(0,  0),
    lrmsd_n_ms  = lrmsd(a1, 0),
    lrmsd_n_ma  = lrmsd(0,  a2),
    lrmsd_n_msa = lrmsd(a1, a2)
  )
}
