# Single-point-mutation scans: generation and preprocessing
# Functions to generate SPM data and reshape it into model-ready arrays

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

#' Generate single-point-mutation (SPM) scan data
#'
#' Runs a single-point-mutation scan over a protein structure: for every site it
#' generates `n_mutations` mutants and records each mutant's measured effects (the
#' energy changes and the per-site / per-mode squared displacement profiles),
#' returning one row per mutant. Mutants are processed one at a time and only the
#' metrics are kept (never a tibble of full mutant structures), so the scan stays
#' memory-efficient on large proteins.
#'
#' @param wt Wild-type protein structure with an elastic network model, as
#'   returned by [setup_enm()].
#' @param n_mutations Number of mutant replicates to generate per site.
#' @param model Name of the mutation model to apply (passed to the ENM machinery).
#' @param sigma Mutation strength (the perturbation magnitude).
#' @param min_sd Minimum sequence separation between coupled sites.
#' @param pdb_site_active Optional integer vector of active-site residue numbers
#'   (PDB numbering).
#' @param seed Optional random seed, for a reproducible scan.
#' @return A tibble with one row per mutant `(j, m)` (`m = 0` is the wild type),
#'   carrying all measured effects of that mutation: scalar energy changes
#'   (`ddg_dv_jm`, `ddg_tds_jm`, `ddgact_dv_jm`, `ddgact_tds_jm`) and list-columns
#'   `site` / `pdb_site` (the site index to PDB-residue map), `dr` (displacement
#'   vector), `dr2_ijm` (per-site squared displacement; each cell a `dr2_i` vector
#'   over response sites `i` for that mutant), `mode` (the normal-mode index
#'   `1:nmodes`), and `dr2_njm` (per-mode squared contribution; each cell a `dr2_n`
#'   vector over response modes `n`).
#' @seealso [setup_enm()] (builds the `wt` input); [preprocess_spm()] and
#'   [preprocess_spm_mode()] (reshape this output for the model).
#' @family spm
#' @examples
#' \dontrun{
#' spm <- generate_spm_data(znb_wt, n_mutations = 10, seed = 1024)
#' head(spm[c("j", "m")])
#' }
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
      # Per-mutant site vector (dr2_i): call penm's primitive directly.
      dr2_i <- penm::delta_structure_dr2i(wt, mut)

      # Mode-form structural divergence: per-normal-mode contribution to dr.
      # Same per-mutant call shape as dr2_i; another reduction of the same mutant.
      dr2_n <- penm::delta_structure_dr2n(wt, mut)

      # Create row and add to results. The list-columns dr2_ijm / dr2_njm carry,
      # per mutant (j,m) row, the (dr2_i) / (dr2_n) vectors.
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
        dr2_ijm = list(dr2_i),
        mode = list(seq_along(dr2_n)),
        dr2_njm = list(dr2_n)
      )

      results <- bind_rows(results, row)

      # Free memory by removing the mutant
      rm(mut)
      if (j %% 50 == 0 && m == n_mutations) gc()  # Occasional garbage collection
    }
  }

  return(results)
}

#' Reshape a mutation scan into model-ready arrays (site form)
#'
#' Converts the raw single-point-mutation scan into the compact arrays the
#' site-form prediction and fitting functions expect. It combines the energy
#' columns into per-mutant stability and activity changes, stacks the per-site
#' squared displacements into a mutant-by-site matrix, and builds the map between
#' the internal site index and the PDB residue number. Wild-type rows are dropped.
#' Run this once and reuse the result across many `(a1, a2)` evaluations.
#'
#' @param spm A single-point-mutation scan, as returned by [generate_spm_data()]
#'   (with the energy columns and the `site`, `pdb_site`, and `dr2_ijm`
#'   list-columns).
#' @return A list with three elements: `energy_data` (a tibble of per-mutant
#'   stability and activity energy changes), `dr2_ijm` (the mutant-by-site matrix
#'   of per-site squared displacements; columns are sites), and `site_map` (a
#'   tibble mapping the site index `i` to its PDB residue number `pdb_site`).
#' @seealso [generate_spm_data()] (produces the input); [calculate_dr2_i_msa()]
#'   and [fit_lrmsd_i_msa_ml()] (consume the output); [preprocess_spm_mode()] for
#'   the mode form.
#' @family spm
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' dim(pp$dr2_ijm)
#' }
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

  # Construct the dr2_ijm matrix incrementally to avoid stack issues
  n_rows <- nrow(spm_filtered)
  n_cols <- length(spm_filtered$site[[1]])
  dr2_ijm <- matrix(0, nrow = n_rows, ncol = n_cols)
  for (i in 1:n_rows) {
    dr2_ijm[i, ] <- spm_filtered$dr2_ijm[[i]]
  }
  colnames(dr2_ijm) <- spm_filtered$site[[1]]  # Site numbers as column names

  # Map internal site index i (= dr2_ijm column names) to structure-anchored
  # pdb_site. The site / pdb_site list-columns are per-row parallel vectors;
  # the first row carries the full ordered mapping.
  site_map <- tibble(
    i = as.integer(spm_filtered$site[[1]]),
    pdb_site = as.integer(spm_filtered$pdb_site[[1]])
  )

  list(energy_data = energy_data, dr2_ijm = dr2_ijm, site_map = site_map)
}

#' Reshape a mutation scan into model-ready arrays (mode form)
#'
#' Mode-indexed counterpart of [preprocess_spm()]: it produces the arrays the
#' mode-form prediction and fitting functions expect. The per-mutant energy table
#' is the same; the difference is that the squared displacements are stacked into a
#' mutant-by-mode matrix rather than mutant-by-site. Because modes are not anchored
#' to residues, there is no PDB-residue map.
#'
#' @param spm A single-point-mutation scan, as returned by [generate_spm_data()]
#'   (with the `mode` and `dr2_njm` list-columns).
#' @return A list with two elements: `energy_data` (a tibble with `j`, `m`,
#'   `ddg_jm`, `ddgact_jm`) and `dr2_njm` (the mutant-by-mode matrix of per-mode
#'   squared displacements; columns are mode indices).
#' @seealso [generate_spm_data()] (produces the input); [calculate_dr2_n_msa()]
#'   and [fit_lrmsd_n_msa_ml()] (consume the output); [preprocess_spm()] for the
#'   site form.
#' @family spm
#' @examples
#' \dontrun{
#' pp <- preprocess_spm_mode(znb_spm)
#' dim(pp$dr2_njm)
#' }
#' @export
preprocess_spm_mode <- function(spm) {
  # Filter out no-mutation cases (m = 0) to align outputs
  spm_filtered <- spm %>% filter(m > 0)

  # Extract energy data (identical to preprocess_spm: weights are axis-agnostic)
  energy_data <- spm_filtered %>%
    mutate(
      ddg_jm = ddg_dv_jm + ddg_tds_jm,
      ddgact_jm = ddgact_dv_jm + ddgact_tds_jm
    ) %>%
    dplyr::select(j, m, ddg_jm, ddgact_jm)

  # Construct the dr2_njm matrix incrementally to avoid stack issues
  n_rows <- nrow(spm_filtered)
  n_cols <- length(spm_filtered$mode[[1]])
  dr2_njm <- matrix(0, nrow = n_rows, ncol = n_cols)
  for (k in 1:n_rows) {
    dr2_njm[k, ] <- spm_filtered$dr2_njm[[k]]
  }
  colnames(dr2_njm) <- spm_filtered$mode[[1]]  # Mode numbers as column names

  list(energy_data = energy_data, dr2_njm = dr2_njm)
}

#' SPM-ensemble averaging weights over mutants
#'
#' Turns the per-mutant fixation probabilities of the MSA model ([pfix_msa()]) into
#' the averaging weights of the single-point-mutation (SPM) ensemble: `weights_jm =
#' pfix_jm / sum(pfix_jm)`, one weight per mutant `(j, m)`, summing to one. The
#' site- and mode-form forward maps ([calculate_dr2_i_msa()], [calculate_dr2_n_msa()])
#' both average their per-mutant divergences with exactly these weights.
#'
#' The normalisation is a property of *this ensemble*, not of the model: it is the
#' SPM scan's estimator of the mutant distribution. A different ensemble (for
#' example one generated by simulating a star tree over several mutations) would turn
#' the same [pfix_msa()] probabilities into different weights. Keeping this step in
#' its own function marks that ensemble-specific seam explicitly.
#'
#' The weights depend on `(a1, a2)` and so must be recomputed whenever the selection
#' strengths change; they cannot be precomputed once from the scan.
#'
#' @param spm_pp Preprocessed single-point-mutation data, the output of
#'   [preprocess_spm()] or [preprocess_spm_mode()] (only its `energy_data` is used,
#'   so either form works).
#' @param a1 Stability selection strength (non-negative).
#' @param a2 Activity selection strength (non-negative).
#' @return A numeric vector of averaging weights over mutants, one per row of
#'   `spm_pp$energy_data`, summing to one.
#' @seealso [pfix_msa()] (the unnormalised model probabilities);
#'   [calculate_dr2_i_msa()] and [calculate_dr2_n_msa()] (average with these weights).
#' @family spm
#' @examples
#' \dontrun{
#' pp <- preprocess_spm(znb_spm)
#' w <- weights_jm_spm(pp, a1 = 1, a2 = 1)
#' sum(w)  # 1
#' }
#' @export
weights_jm_spm <- function(spm_pp, a1, a2) {
  energy_data <- spm_pp$energy_data
  pfix_jm <- pfix_msa(energy_data$ddg_jm, energy_data$ddgact_jm, a1, a2)
  pfix_jm / sum(pfix_jm)
}
