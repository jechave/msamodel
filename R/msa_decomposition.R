#' Split a divergence profile into mutation, stability, and activity contributions
#'
#' Decomposes a structural-divergence profile into the part attributable to
#' mutation alone and the parts added by each selection pressure. The split is
#' sequential along the model progression no-selection -> stability -> activity:
#' `phi_mut` is the mutation-only divergence (the MM variant), `phi_stab` is the
#' extra divergence explained when stability selection is added (MS minus MM), and
#' `phi_act` is the further change when activity selection is added (MSA minus MS).
#' The three contributions sum exactly to the full-model profile.
#'
#' The inputs are four numeric vectors of equal length -- the divergence under each
#' of the four model variants. The function works on the values alone and does not
#' care whether they are indexed by site or by mode, so the same call serves both;
#' the caller supplies whichever four columns it holds and attaches the result.
#'
#' The public entry points that own the four-variant evaluation and attach the
#' result are [calculate_decomposition_i_msa()] / [calculate_decomposition_n_msa()]
#' (forward, at one `(a1, a2)`) and [calculate_decomposition_samples()] (across
#' posterior samples); this kernel is the pure math they share and is internal.
#'
#' @param mm,ms,ma,msa Numeric vectors of equal length giving the divergence under
#'   the four model variants: mutation-only (`mm`), plus-stability (`ms`),
#'   plus-activity (`ma`), and the full model (`msa`). `ma` is accepted but not used
#'   by the current sequential formula; it is kept so an alternative decomposition
#'   that needs it can be added later without changing the signature.
#' @return A named list of three numeric vectors -- `phi_mut`, `phi_stab`,
#'   `phi_act` -- one value per input element.
#' @noRd
calculate_msa_decomposition <- function(mm, ms, ma, msa) {
  # Sequential M0 -> MM -> MS -> MSA decomposition. `ma` is unused here but kept
  # in the signature so the future Shapley `method` is a non-breaking addition.
  list(
    phi_mut  = mm,
    phi_stab = ms - mm,
    phi_act  = msa - ms
  )
}

#' Apply the divergence decomposition to every posterior sample
#'
#' Runs `calculate_msa_decomposition()` on each row of a table of per-sample,
#' per-site divergence profiles, turning the four model-variant columns into the
#' three contribution columns (`phi_mut`, `phi_stab`, `phi_act`). The split is
#' computed independently for each row, so samples and sites need no grouping.
#'
#' @param data A data frame with columns `sample_id`, `i`, and the four model
#'   variants `lrmsd_i_mm`, `lrmsd_i_ms`, `lrmsd_i_ma`, `lrmsd_i_msa`. A `pdb_site`
#'   column, if present, is carried through.
#' @return A data frame with `sample_id`, `i` (and `pdb_site` if it was present),
#'   and the three contribution columns `phi_mut`, `phi_stab`, `phi_act`.
#' @seealso `calculate_msa_decomposition()` (the per-row split it applies);
#'   [calculate_decomposition_summary()] (summarises the result across samples).
#' @family decomposition
#' @examples
#' \dontrun{
#' analysis <- run_msa_bayesian_analysis(znb_spm, znb_profile)
#' phi <- calculate_decomposition_samples(analysis$prediction_samples)
#' head(phi)
#' }
#' @export
calculate_decomposition_samples <- function(data) {
  # Validate input columns
  required_cols <- c("sample_id", "i", "lrmsd_i_mm", "lrmsd_i_ms", "lrmsd_i_ma", "lrmsd_i_msa")

  missing_cols <- required_cols[!required_cols %in% colnames(data)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  keep <- c("sample_id", "i", if ("pdb_site" %in% colnames(data)) "pdb_site")
  phi <- calculate_msa_decomposition(
    data$lrmsd_i_mm, data$lrmsd_i_ms, data$lrmsd_i_ma, data$lrmsd_i_msa
  )
  result <- dplyr::bind_cols(data[keep], tibble::as_tibble(phi))

  return(result)
}

#' Summarise the divergence decomposition across posterior samples
#'
#' Reduces the per-sample contributions to a posterior summary, with one row per
#' site and contribution. For each contribution (`phi_mut`, `phi_stab`,
#' `phi_act`) it reports the mean, standard deviation, median, and 95% credible
#' interval over the samples. A `pdb_site` column is carried through when present.
#'
#' @param decomposition_samples A data frame with columns `sample_id`, `i`, and the
#'   three contribution columns `phi_mut`, `phi_stab`, `phi_act` (as returned by
#'   [calculate_decomposition_samples()]); `pdb_site` is carried through if present.
#' @return A long-format data frame with `i` (and `pdb_site` if present), a
#'   `component` column naming the contribution, and the summary columns `mean`,
#'   `sd`, `median`, `lower`, and `upper`.
#' @seealso [calculate_decomposition_samples()] (produces the input);
#'   [run_msa_bayesian_analysis()] (which calls this as its final step).
#' @family decomposition
#' @examples
#' \dontrun{
#' analysis <- run_msa_bayesian_analysis(znb_spm, znb_profile)
#' summ <- calculate_decomposition_summary(analysis$decomposition_samples)
#' head(summ)
#' }
#' @export
calculate_decomposition_summary <- function(decomposition_samples) {
  # Validate input columns
  required_cols <- c("sample_id", "i", "phi_mut", "phi_stab", "phi_act")

  missing_cols <- required_cols[!required_cols %in% colnames(decomposition_samples)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Group by the site key(s): carry pdb_site through when present
  group_keys <- if ("pdb_site" %in% colnames(decomposition_samples)) {
    c("i", "pdb_site", "component")
  } else {
    c("i", "component")
  }

  # Calculate summary statistics in long format
  result <- decomposition_samples %>%
    # Step 1: Convert to long format first
    pivot_longer(
      cols = starts_with("phi_"),
      names_to = "component",
      values_to = "value"
    ) %>%
    # Step 2: Group by residue and component
    group_by(across(all_of(group_keys))) %>%
    # Step 3: Calculate summary statistics
    summarise(
      mean = mean(value, na.rm = TRUE),
      sd = sd(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      lower = quantile(value, 0.025, na.rm = TRUE),
      upper = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  return(result)
}
