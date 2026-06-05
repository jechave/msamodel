#' MSA Data Preparation
#' Functions to preprocess SPM data for Mutation Stability Activity (MSA) analysis

#' Preprocess SPM data for matrix-based Bayesian analysis
#'
#' @param spm SPM tibble with columns j, m, ddg_dv_jm, ddg_tds_jm, ddgact_dv_jm, ddgact_tds_jm, and list-columns site and dr2
#' @return List containing energy_data and dr2mat
#' @family spm
#' @export
preprocess_spm <- function(spm) {
  # Filter out no-mutation cases (m = 0) to align outputs
  spm_filtered <- spm %>% filter(m > 0)

  # Extract energy data
  energy_data <- spm_filtered %>%
    mutate(
      ddg_jm = ddg_dv_jm + ddg_tds_jm,
      ddgact_jm = ddgact_dv_jm + ddgact_tds_jm
    ) %>%
    dplyr::select(j, m, ddg_jm, ddgact_jm)

  # Construct dr2 matrix incrementally to avoid stack issues
  n_rows <- nrow(spm_filtered)
  n_cols <- length(spm_filtered$site[[1]])
  dr2mat <- matrix(0, nrow = n_rows, ncol = n_cols)
  for (i in 1:n_rows) {
    dr2mat[i, ] <- spm_filtered$dr2[[i]]
  }
  colnames(dr2mat) <- spm_filtered$site[[1]]  # Site numbers as column names

  list(energy_data = energy_data, dr2mat = dr2mat)
}
