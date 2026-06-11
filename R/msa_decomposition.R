#' Calculate the phi (site-level) decomposition of the MSA model for the lrmsd metric
#'
#' Splits the per-site log structural divergence into mutation, stability-
#' selection, and activity-selection components (phi_mut, phi_stab, phi_act).
#' (Earlier versions called these `shap_*` / a "Shapley" decomposition; that was
#' a misnomer and the columns were renamed to `phi_*` in v0.2.)
#'
#' @param data Data frame with columns i, lrmsd_mm, lrmsd_ms, lrmsd_ma, lrmsd_msa
#'   (and optionally pdb_site, which is carried through if present)
#' @return Data frame with the site key(s) (i, and pdb_site if present) and the
#'   phi decomposition columns (phi_mut, phi_stab, phi_act)
#' @family decomposition
#' @export
calculate_msa_decomposition <- function(data) {
  # Validate input columns
  required_cols <- c("i", "lrmsd_mm", "lrmsd_ms", "lrmsd_ma", "lrmsd_msa")

  missing_cols <- required_cols[!required_cols %in% colnames(data)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Calculate phi components (carry pdb_site through when present)
  result <- data %>%
    transmute(
      i = i,
      phi_mut = lrmsd_mm,
      phi_stab = 0.5 * (lrmsd_ms - lrmsd_mm + lrmsd_msa - lrmsd_ma),
      phi_act = 0.5 * (lrmsd_ma - lrmsd_mm + lrmsd_msa - lrmsd_ms)
    )
  if ("pdb_site" %in% colnames(data)) {
    result <- dplyr::bind_cols(result["i"], data["pdb_site"], result[c("phi_mut", "phi_stab", "phi_act")])
  }

  return(result)
}

#' Calculate the MSA phi decomposition for multiple samples
#'
#' @param data Data frame with columns sample_id, i, lrmsd_mm, lrmsd_ms, lrmsd_ma, lrmsd_msa
#' @return Data frame with only sample_id, i, and the phi decomposition columns
#' @family decomposition
#' @export
calculate_decomposition_samples <- function(data) {
  # Validate input columns
  required_cols <- c("sample_id", "i", "lrmsd_mm", "lrmsd_ms", "lrmsd_ma", "lrmsd_msa")

  missing_cols <- required_cols[!required_cols %in% colnames(data)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Process each sample group
  result <- data %>%
    group_by(sample_id) %>%
    group_split() %>%
    purrr::map_dfr(function(sample_df) {
      # Calculate decomposition for this sample
      decomp_values <- calculate_msa_decomposition(sample_df)

      # Add sample_id back for identification
      decomp_values %>%
        mutate(sample_id = unique(sample_df$sample_id))
    })

  # Ensure columns are properly ordered
  result <- result %>%
    select(sample_id, i, everything())

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
