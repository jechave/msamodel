#' Get displacement vector dr between two structures
#'
#' @param wt Wild type protein structure
#' @param mut Mutant protein structure
#' @return Displacement vector dr
#' @noRd
delta_structure_dr <- function(wt, mut) {
  # Check that sites match
  stopifnot(all(wt$node$pdb_site == mut$node$pdb_site)) # no indels
  stopifnot(all(wt$node$site == mut$node$site)) # no indels

  # Return the raw displacement vector
  mut$node$xyz - wt$node$xyz
}

#' Calculate squared displacement per site
#'
#' @param wt Wild type protein structure
#' @param mut Mutant protein structure
#' @return Vector of squared displacements per site
#' @noRd
delta_structure_dr2 <- function(wt, mut) {
  # Check that sites match
  stopifnot(all(wt$node$pdb_site == mut$node$pdb_site)) # no indels
  stopifnot(all(wt$node$site == mut$node$site)) # no indels

  # Calculate displacement vector
  dr <- mut$node$xyz - wt$node$xyz

  # Calculate squared displacement per site
  n_sites <- length(wt$node$site)
  dr2 <- numeric(n_sites)

  for (i in 1:n_sites) {
    # Get indices for the x, y, z coordinates of site i
    idx <- (3*i-2):(3*i)
    # Calculate squared displacement
    dr2[i] <- sum(dr[idx]^2)
  }

  return(dr2)
}

#' Process protein mutants one by one and calculate only necessary metrics
#'
#' A memory-efficient approach to generating and analyzing protein mutations.
#' Instead of creating a large tibble with all mutant structures, this function
#' processes each mutation individually, calculates the required metrics, and
#' only stores the results.
#'
#' @param wt Wild type protein structure with ENM
#' @param n_mutations Number of mutations to generate per site
#' @param model Mutation model to use
#' @param sigma Mutation strength parameter
#' @param min_sd Minimum sequence distance
#' @param pdb_site_active Active site residue numbers (optional)
#' @param seed Random seed for reproducibility
#' @return Tibble with one row per mutant `(j, m)` (`m = 0` is the wild type),
#'   carrying all measured effects of that mutation: scalar energy changes
#'   (`ddg_dv_jm`, `ddg_tds_jm`, `ddgact_dv_jm`, `ddgact_tds_jm`) and list-columns
#'   `site` / `pdb_site` (the site index <-> PDB-residue map), `dr` (displacement
#'   vector), `dr2` (per-site squared displacement), `mode` (the normal-mode index
#'   `1:nmodes`), and `dr2n` (per-mode squared contribution to `dr`).
#' @family spm
#' @export
generate_spm_data <- function(wt, n_mutations = 10,
                             model = "lfenm", sigma = 0.3,
                             min_sd = 2, pdb_site_active = NULL,
                             seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Get site information once
  site_vector <- get_site(wt)
  pdb_site_vector <- get_pdb_site(wt)
  n_sites <- length(site_vector)

  # Initialize empty results tibble
  results <- tibble()

  # Process each site and mutation state
  for (j in site_vector) {
    # Create a progress indicator
    if (j %% 10 == 0) cat("Processing site", j, "of", n_sites, "\n")

    for (m in 0:n_mutations) {
      # Generate single mutant (m=0 is wild type)
      mut <- get_mutant_site(wt, j, m,
                             mut_model = model,
                             mut_dl_sigma = sigma,
                             mut_sd_min = min_sd,
                             seed = seed)

      # Calculate energies
      ddg_dv_jm <- ddg_dv(wt, mut)
      ddg_tds_jm <- ddg_tds(wt, mut)

      # Calculate activity energies if requested
      if (!is.null(pdb_site_active)) {
        ddgact_dv_jm <- ddgact_dv(wt, mut, pdb_site_active = pdb_site_active)
        ddgact_tds_jm <- ddgact_tds(wt, mut, pdb_site_active = pdb_site_active)
      } else {
        ddgact_dv_jm <- NA_real_
        ddgact_tds_jm <- NA_real_
      }

      # Calculate displacement vectors
      dr <- delta_structure_dr(wt, mut)
      dr2 <- delta_structure_dr2(wt, mut)

      # Mode-form structural divergence: per-normal-mode contribution to dr.
      # Same per-mutant call shape as dr2; another reduction of the same mutant.
      dr2n <- penm::delta_structure_dr2n(wt, mut)

      # Create row and add to results
      row <- tibble(
        j = j,
        m = m,
        ddg_dv_jm = ddg_dv_jm,
        ddg_tds_jm = ddg_tds_jm,
        ddgact_dv_jm = ddgact_dv_jm,
        ddgact_tds_jm = ddgact_tds_jm,
        site = list(site_vector),
        pdb_site = list(pdb_site_vector),
        dr = list(dr),
        dr2 = list(dr2),
        mode = list(seq_along(dr2n)),
        dr2n = list(dr2n)
      )

      results <- bind_rows(results, row)

      # Free memory by removing the mutant
      rm(mut)
      if (j %% 50 == 0 && m == n_mutations) gc()  # Occasional garbage collection
    }
  }

  return(results)
}
