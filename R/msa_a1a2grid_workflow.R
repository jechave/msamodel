#' MSA A1A2 Grid Workflow
#' Site-specific structural analysis over a1, a2 grid

#' Aggregate site-specific structural differences with MSA weighting
#'
#' @param spm SPM tibble
#' @param a1a2grid Tibble with a1, a2
#' @param verbose Logical
#' @return Tibble with site-specific results
#' @family grid
#' @export
calculate_dr2i_msa_a1a2grid <- function(spm, a1a2grid, verbose = TRUE) {
  if (verbose) cat("Preparing site-specific data from SPM tibble\n")
  spm_energies_and_dr2mat <- preprocess_spm(spm)
  site_ids <- as.integer(colnames(spm_energies_and_dr2mat$dr2mat))

  if (verbose) cat("Calculating for", nrow(a1a2grid), "selection parameter combinations\n")
  result <- map_dfr(1:nrow(a1a2grid), function(idx) {
    if (verbose && idx %% 5 == 0)
      cat("  Processing parameter set", idx, "of", nrow(a1a2grid), "\n")
    a1 <- a1a2grid$a1[idx]
    a2 <- a1a2grid$a2[idx]

    # Combined selection
    dr2_msa <- calculate_dr2i_msa(spm_energies_and_dr2mat, a1, a2)$dr2_i

    # Stability only
    dr2_ms <- calculate_dr2i_msa(spm_energies_and_dr2mat, a1, 0)$dr2_i

    # Activity only
    dr2_ma<- calculate_dr2i_msa(spm_energies_and_dr2mat, 0, a2)$dr2_i

    # No selection
    dr2_mm <- colMeans(spm_energies_and_dr2mat$dr2mat)

    tibble(
      a1 = a1,
      a2 = a2,
      i = site_ids,
      dr2_msa,
      dr2_mm,
      dr2_ms,
      dr2_ma
    )
  })

  return(result)
}

#' Define Grid of Selection Parameters for MSA Analysis
#'
#' Creates a tibble with combinations of stability (a1) and activity (a2) selection parameters
#' for systematic exploration of parameter space in the Mutation-Stability-Activity model.
#'
#' @param a1_max Maximum value for stability selection parameter (minimum is always 0)
#' @param log2_a2_plus1_max Maximum value for log2(a2+1) transformation
#' @param n_a1 Number of grid points for a1 dimension (including 0)
#' @param n_a2 Number of grid points for a2 dimension (including 0)
#' @return Tibble with columns a1 and a2, containing parameter combinations
#' @family grid
#' @export
define_selection_grid <- function(a1_max = 10,
                                  log2_a2_plus1_max = 12,
                                  n_a1 = 6,
                                  n_a2 = 9) {

  # Generate a1 values with linear spacing
  if (n_a1 <= 2) {
    a1_values <- c(0, a1_max)
  } else {
    a1_values <- seq(0, a1_max, length.out = n_a1)
  }

  # Generate log2(a2+1) values, then transform to a2 values
  if (n_a2 <= 2) {
    log2_a2_plus1_values <- c(0, log2_a2_plus1_max)
  } else {
    # Linearly spaced in log2(a2+1) space
    log2_a2_plus1_values <- seq(0, log2_a2_plus1_max, length.out = n_a2)
  }

  # Convert log2(a2+1) values to a2 values
  a2_values <- 2^log2_a2_plus1_values - 1

  # Create all combinations
  grid <- expand.grid(a1 = a1_values, a2 = a2_values) %>%
    as_tibble()

  # Round values to make them more readable
  grid <- grid %>%
    mutate(
      a1 = round(a1, 4),
      a2 = round(a2, 4)
    )

  # Sort the grid for consistent results
  grid <- grid %>% arrange(a1, a2)

  return(grid)
}
