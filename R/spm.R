# Single-point-mutation scans: generation and preprocessing
# Functions to generate SPM data and reshape it into model-ready arrays

#' Run the raw single-point-mutation scan (per-mutant record)
#'
#' The lossless SPM primitive: for every site it generates `n_mutations` mutants and
#' records each mutant's measured effects, returning one **row per mutant**. Mutants are
#' processed one at a time and only the metrics are kept (never a tibble of full mutant
#' structures), so the scan stays memory-efficient on large proteins. This is the
#' internal core that [generate_spm()] wraps: it keeps the per-mutant record of
#' reductions, which [generate_spm()] then reshapes into the lean, model-ready `spm`
#' object. A future motion arm adds its quantities here the same way `dr2_i` / `dr2_n`
#' are added -- by calling penm's own per-mutant primitives (`penm::delta_motion_*`,
#' which take `(wt, mut)`) inside this loop, while the mutant is live.
#'
#' The list-column names here keep the index-signature convention (`dr2_ijm`,
#' `dr2_njm`): this is internal, and the letters state the shape exactly -- one row per
#' mutant `(j, m)`, each cell a vector over response sites `i` (or modes `n`). The
#' public object renames them to `dr2mat_site` / `dr2mat_mode`; [generate_spm()]
#' is the boundary between the two vocabularies.
#'
#' @param wt Wild-type protein structure with an elastic network model, as
#'   returned by [set_enm()].
#' @param n_mutations Number of mutant replicates to generate per site.
#' @param model Name of the mutation model to apply (passed to the ENM machinery).
#' @param sigma Mutation strength (the perturbation magnitude).
#' @param min_sd Minimum sequence separation between coupled sites.
#' @param pdb_site_active Optional integer vector of active-site residue numbers
#'   (PDB numbering).
#' @param ensemble Which sample of mutants to draw; see [generate_spm()].
#' @return A tibble with one row per mutant `(j, m)` (`m = 0` is the wild type),
#'   carrying all measured effects of that mutation: scalar energy changes
#'   (`ddg_dv_jm`, `ddg_tds_jm`, `ddgact_dv_jm`, `ddgact_tds_jm`) and list-columns
#'   `site` / `pdb_site` (the site index to PDB-residue map), `dr2_ijm` (per-site
#'   squared displacement; each cell a `dr2_i` vector over response sites `i` for that
#'   mutant), `mode` (penm's normal-mode index, [penm::get_mode()]), and `dr2_njm`
#'   (per-mode squared contribution; each cell a `dr2_n` vector over response modes `n`).
#' @noRd
generate_spm_core <- function(wt, n_mutations = 10,
                             model = "lfenm", sigma = 0.3,
                             min_sd = 2, pdb_site_active = NULL,
                             ensemble) {
  # Get site information once
  site_vector <- get_site(wt)
  pdb_site_vector <- get_pdb_site(wt)
  # The mode index is penm's own (`nma$mode`), not a locally counted one: asking penm
  # means the label cannot drift from the modes it names if penm's mode set ever changes.
  mode_vector <- penm::get_mode(wt)

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
                             ensemble = ensemble)

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

      # Two reductions of the same mutant displacement, both from penm's primitives:
      # per-site (dr2_i) and per-normal-mode (dr2_n).
      dr2_i <- penm::delta_structure_dr2i(wt, mut)
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
        dr2_ijm = list(dr2_i),
        mode = list(mode_vector),
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

#' Run a single-point-mutation scan over a protein
#'
#' Runs a single-point-mutation (SPM) scan: every site of the protein is mutated
#' `n_mutations` times, and for each mutant the scan records two free-energy changes
#' and how far every residue moved.
#'
#' Mutating site `j` perturbs the equilibrium length of every network edge touching
#' `j`, each by an independent draw from `Normal(0, sigma)` (only edges whose sequence
#' separation is at least `min_sd` are perturbed). Those altered lengths put the
#' structure under strain, and relaxing it displaces the residues. The scan stores the
#' squared displacements.
#'
#' @section The mutation-response matrices:
#'
#' `dr2mat_site` and `dr2mat_mode` both have **one row per mutant**, in scan order;
#' row `k` of either is the mutant described by row `k` of `energy_data`. They differ
#' in what a column is:
#'
#' - `dr2mat_site` — one column per **residue**. Cell `[k, i]` is the squared
#'   displacement of residue `i` (Å²) under mutant `k`.
#' - `dr2mat_mode` — one column per **normal mode**. Cell `[k, n]` is the squared
#'   projection of the same displacement onto mode `n`.
#'
#' The modes are orthonormal, so each row sums to the same total in both matrices: they
#' hold one displacement written in two bases. Neither carries column names — the
#' response index is the column position — so that `colSums()` and the like return bare
#' vectors.
#'
#' @param wt The protein to mutate, carrying its elastic network model.
#' @param n_mutations Mutants generated per site. The scan produces
#'   `n_sites * n_mutations` mutants in all, which is the row count of both matrices.
#' @param model Mutation model, passed to penm as `mut_model`: `"lfenm"` (the default)
#'   or `"sclfenm"`. Any other value is an error.
#' @param sigma Standard deviation (Å) of the `Normal(0, sigma)` perturbation applied
#'   to each affected edge length. It sets the size of a mutation.
#' @param min_sd Minimum sequence separation for an edge to be perturbed: at the
#'   default `2`, edges to immediate sequence neighbours are left alone.
#' @param pdb_site_active Active-site residue numbers, in **PDB numbering**. If `NULL`,
#'   no activity energy is computed and `energy_data$ddgact` is `NA` throughout, which
#'   leaves the model's activity term unusable.
#' @param ensemble Which sample of mutants to draw. The mutations are random, so a scan
#'   is one of many possible ones; a given `ensemble` always yields the same mutants, on
#'   any machine, and a different one yields an equally valid different set. Use a single
#'   value across anything you intend to compare, and record it -- it is what makes the
#'   scan reproducible. Any integer. See [penm::penm_ensemble].
#' @return An `spm` object (a classed list) with five elements:
#'   \describe{
#'     \item{`energy_data`}{One row per mutant, aligned with the matrix rows: `j` (the
#'       mutated site), `m` (which replicate), `ddg` (change in folding free energy)
#'       and `ddgact` (change in activity free energy, or `NA` — see
#'       `pdb_site_active`).}
#'     \item{`dr2mat_site`}{`[mutant x residue]` squared displacements (Å²).}
#'     \item{`dr2mat_mode`}{`[mutant x mode]` squared displacements.}
#'     \item{`site_map`}{The PDB residue number (`pdb_site`) for each `dr2mat_site`
#'       column position (`site`).}
#'     \item{`mode_map`}{The mode index for each `dr2mat_mode` column.}
#'   }
#' @seealso [set_enm()] (builds the elastic network model this needs);
#'   [calculate_profiles()] and [fit_lrmsd_msa_site()] (consume the returned object).
#' @family spm
#' @examples
#' \dontrun{
#' ex  <- function(f) system.file("extdata", f, package = "msamodel")
#' wt  <- set_enm(bio3d::read.pdb(ex("1znb_A.pdb")), node = "ca",
#'                       model = "ming_wall", d_max = 10.5, frustrated = FALSE)
#' act <- readr::read_csv(ex("znb_active_site.csv"))
#' spm <- generate_spm(wt, n_mutations = 10, pdb_site_active = act$pdb_site,
#'                     ensemble = 1L)
#'
#' dim(spm$dr2mat_site)   # mutants x residues
#' dim(spm$dr2mat_mode)   # mutants x modes
#'
#' # The same displacement in two bases: each row sums to the same total.
#' all.equal(rowSums(spm$dr2mat_site), rowSums(spm$dr2mat_mode))
#' }
#' @export
generate_spm <- function(wt, n_mutations = 10,
                             model = "lfenm", sigma = 0.3,
                             min_sd = 2, pdb_site_active = NULL,
                             ensemble) {
  scan <- generate_spm_core(wt, n_mutations = n_mutations, model = model,
                            sigma = sigma, min_sd = min_sd,
                            pdb_site_active = pdb_site_active, ensemble = ensemble)

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
