# MSA model forward map (divergence calculators at given (a1, a2))
# Functions to compute the model's predicted structural-divergence profiles

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
#' Pass the four returned columns to [calculate_msa_decomposition()] to split the
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
#'   [calculate_msa_decomposition()] (splits the four columns into contributions);
#'   [calculate_lrmsd_n_nested_models()] (the mode-indexed counterpart).
#' @family model
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' head(calculate_lrmsd_i_nested_models(pp, a1 = 1, a2 = 1))
#' }
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
#' [calculate_msa_decomposition()] to split the divergence into its contributions.
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
  lrmsd <- function(p1, p2) log(sqrt(calculate_dr2_n_msa(spm_pp, p1, p2)$dr2_n))

  tibble(
    n           = as.integer(colnames(spm_pp$dr2_njm)),
    lrmsd_n_mm  = lrmsd(0,  0),
    lrmsd_n_ms  = lrmsd(a1, 0),
    lrmsd_n_ma  = lrmsd(0,  a2),
    lrmsd_n_msa = lrmsd(a1, a2)
  )
}
