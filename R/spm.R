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

#' Run the raw single-point-mutation scan (per-mutant record)
#'
#' The lossless SPM primitive: for every site it generates `n_mutations` mutants and
#' records each mutant's measured effects, returning one **row per mutant**. Mutants are
#' processed one at a time and only the metrics are kept (never a tibble of full mutant
#' structures), so the scan stays memory-efficient on large proteins. This is the
#' internal core that [generate_spm_data()] wraps: it keeps the full record -- including
#' the raw Cartesian displacement `dr` -- so a future Cartesian/motion calculation can
#' obtain it without re-running the scan. The public [generate_spm_data()] reshapes this
#' record into the lean, model-ready `spm` object and does not surface `dr`.
#'
#' The list-column names here keep the index-signature convention (`dr2_ijm`,
#' `dr2_njm`): this is internal, and the letters state the shape exactly -- one row per
#' mutant `(j, m)`, each cell a vector over response sites `i` (or modes `n`). The
#' public object renames them to `dr2mat_site` / `dr2mat_mode`; [generate_spm_data()]
#' is the boundary between the two vocabularies.
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
#' @noRd
generate_spm_core <- function(wt, n_mutations = 10,
                             model = "lfenm", sigma = 0.3,
                             min_sd = 2, pdb_site_active = NULL,
                             seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Get site information once
  site_vector <- get_site(wt)
  pdb_site_vector <- get_pdb_site(wt)

  # Initialize empty results tibble
  results <- tibble()

  # Process each site and mutation state
  for (j in site_vector) {
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

#' Generate single-point-mutation (SPM) scan data
#'
#' Runs a single-point-mutation scan over a protein structure and returns the
#' **model-ready** `spm` object: for every site it generates `n_mutations` mutants,
#' records each mutant's energy changes and squared structural divergence, and packs the
#' result into the compact arrays every `calculate_*`, `fit_*`, and `predict_*` function
#' consumes. This is the single object the rest of the package works from -- there is no
#' separate preprocessing step.
#'
#' The returned object carries, computed once, both response axes at the same time: the
#' per-site divergence matrix `dr2_ijm` (each site value is the squared displacement of a
#' residue) and the per-mode divergence matrix `dr2_njm` (each value is the squared
#' contribution of a normal mode). A calculation picks the axis it needs; the object does
#' not have to be regenerated to switch between them.
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
#' @return An `spm` object: a list with five elements -- `energy_data` (a tibble of
#'   per-mutant stability and activity energy changes, `j`, `m`, `ddg`, `ddgact`,
#'   wild-type rows dropped), `dr2mat_site` (the mutant-by-site squared-divergence
#'   matrix), `dr2mat_mode` (the mutant-by-mode squared-divergence matrix), `site_map`
#'   (a tibble mapping the internal site index `site` to its PDB residue number
#'   `pdb_site`), and `mode_map` (a tibble carrying the mode index `mode`).
#' @seealso [setup_enm()] (builds the `wt` input); [calculate_profiles()] and
#'   [fit_lrmsd_msa_site()] (consume the returned object).
#' @family spm
#' @examples
#' \dontrun{
#' ex  <- function(f) system.file("extdata", f, package = "msamodel")
#' wt  <- setup_enm(bio3d::read.pdb(ex("1znb_A.pdb")), node = "ca", d_max = 10.5)
#' act <- readr::read_csv(ex("znb_active_site.csv"))
#' spm <- generate_spm_data(wt, n_mutations = 10, pdb_site_active = act$pdb_site,
#'                          seed = 1024)
#' dim(spm$dr2mat_site)
#' dim(spm$dr2mat_mode)
#' }
#' @export
generate_spm_data <- function(wt, n_mutations = 10,
                             model = "lfenm", sigma = 0.3,
                             min_sd = 2, pdb_site_active = NULL,
                             seed = NULL) {
  scan <- generate_spm_core(wt, n_mutations = n_mutations, model = model,
                            sigma = sigma, min_sd = min_sd,
                            pdb_site_active = pdb_site_active, seed = seed)

  # Drop the wild-type rows (m = 0): every model-ready array is per MUTANT.
  scan <- scan %>% filter(m > 0)

  # --- per-mutant energies. Computed once: the fixation weights are axis-agnostic,
  # so the site and mode arrays below share this one table.
  energy_data <- scan %>%
    mutate(
      ddg = ddg_dv_jm + ddg_tds_jm,
      ddgact = ddgact_dv_jm + ddgact_tds_jm
    ) %>%
    dplyr::select(j, m, ddg, ddgact)

  n_rows <- nrow(scan)

  # --- site axis: stack the per-mutant dr2_i vectors into [mutant x site].
  # Built row-at-a-time rather than with do.call(rbind, ...) to avoid stack issues.
  n_sites <- length(scan$site[[1]])
  dr2mat_site <- matrix(0, nrow = n_rows, ncol = n_sites)
  for (i in 1:n_rows) {
    dr2mat_site[i, ] <- scan$dr2_ijm[[i]]
  }

  # --- mode axis: the same stacking, one column per normal mode.
  n_modes <- length(scan$mode[[1]])
  dr2mat_mode <- matrix(0, nrow = n_rows, ncol = n_modes)
  for (k in 1:n_rows) {
    dr2mat_mode[k, ] <- scan$dr2_njm[[k]]
  }

  # Map the internal site index to the structure-anchored pdb_site. The index is the
  # dr2mat_site column position (site[[1]] is 1:n_sites by construction); it is carried
  # here explicitly rather than as matrix colnames (which would leak onto colSums
  # results). The site / pdb_site list-columns are per-row parallel vectors; the first
  # row carries the full ordered mapping.
  #
  # The column is `site`, not `i`: this map translates EITHER a response site (`i`) or a
  # mutated site (`j`) to its PDB label -- it is indifferent to the role. `i`/`j` keep
  # their role-specific meaning in the scan arrays, where the contrast is real.
  site_map <- tibble(
    site = as.integer(scan$site[[1]]),
    pdb_site = as.integer(scan$pdb_site[[1]])
  )

  # Map the mode index to its dr2mat_mode column position. The mode-axis counterpart of
  # site_map, and built the same way: taken FROM THE SCAN (mode[[1]], the per-row
  # parallel vector whose first row carries the full ordered index), not manufactured
  # here. Carried explicitly as a one-column tibble rather than as matrix colnames,
  # which would leak onto colSums results. Modes are not residue-anchored, so there is
  # no pdb_ analogue -- the index is the whole map.
  mode_map <- tibble(mode = as.integer(scan$mode[[1]]))

  spm <- list(
    energy_data = energy_data,
    dr2mat_site = dr2mat_site,
    dr2mat_mode = dr2mat_mode,
    site_map    = site_map,
    mode_map    = mode_map
  )
  class(spm) <- "spm"
  spm
}
