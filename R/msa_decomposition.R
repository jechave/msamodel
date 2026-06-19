#' Phi decomposition of an MSA lrmsd profile (pure, axis-agnostic)
#'
#' Splits a log structural-divergence profile into mutation, stability-selection,
#' and activity-selection components (phi_mut, phi_stab, phi_act). This is the
#' **sequential / additive** decomposition along the model progression
#' M0 -> MM -> MS -> MSA: each term is the incremental divergence added by turning
#' on the next selection pressure. The terms telescope, so
#' phi_mut + phi_stab + phi_act = msa.
#'
#' This is pure math on four vectors: it takes the four nested-model lrmsd values
#' and returns the three phi vectors. It knows nothing about sites vs modes,
#' column names, or keys, so it serves both the site (`i`) and mode (`n`) axes
#' unchanged -- the caller supplies whichever four columns it holds (e.g.
#' `lrmsd_i_*` or, later, `lrmsd_n_*`) and attaches the result.
#'
#' (A different, symmetric "Shapley" decomposition of the same four variants also
#' exists; an earlier version computed that one by mistake. The sequential form
#' here is the one used for the paper analysis.)
#'
#' @param mm,ms,ma,msa Numeric vectors: the four nested-model lrmsd values
#'   (mutation-only, +stability, +activity, full). `ma` is not used by the
#'   sequential formula but is kept as an argument for forward-compatibility with
#'   a future `method` switch to the Shapley variant (which needs it).
#' @return A named list of three numeric vectors: `phi_mut`, `phi_stab`,
#'   `phi_act`. Splice into a tibble with
#'   `dplyr::mutate(df, !!!calculate_msa_decomposition(mm, ms, ma, msa))`.
#' @family decomposition
#' @export
calculate_msa_decomposition <- function(mm, ms, ma, msa) {
  # Sequential M0 -> MM -> MS -> MSA decomposition. `ma` is unused here but kept
  # in the signature so the future Shapley `method` is a non-breaking addition.
  list(
    phi_mut  = mm,
    phi_stab = ms - mm,
    phi_act  = msa - ms
  )
}

#' Calculate the MSA phi decomposition for multiple samples
#'
#' Per-row (per sample x site) sequential phi decomposition, built by splicing the
#' pure [calculate_msa_decomposition()] over the four nested-model columns. The
#' formula is row-wise, so no per-sample grouping is needed.
#'
#' @param data Data frame with columns sample_id, i, lrmsd_i_mm, lrmsd_i_ms,
#'   lrmsd_i_ma, lrmsd_i_msa (and optionally pdb_site, carried through if present)
#' @return Data frame with sample_id, i (and pdb_site if present), and the phi
#'   decomposition columns (phi_mut, phi_stab, phi_act)
#' @family decomposition
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

#' Summarize the MSA phi decomposition across samples in long format
#'
#' @param decomposition_samples Data frame with sample_id, i, phi_mut, phi_stab,
#'   phi_act columns (and optionally pdb_site, carried through if present)
#' @return Data frame with i (and pdb_site if present), component, and summary
#'   statistics in long format
#' @family decomposition
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
