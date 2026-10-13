# Single-point-mutation scans: generation and preprocessing
# Functions to generate SPM data and reshape it into model-ready arrays

#' Run the raw single-point-mutation scan (per-mutant record)
#'
#' The lossless SPM primitive: it gives every site each letter of the alphabet other
#' than the wild type's, one mutant per letter, and records each mutant's measured
#' effects, returning one **row per mutant**. Mutants are processed one at a time and
#' only the metrics are kept (never a tibble of full mutant structures), so the scan
#' stays memory-efficient on large proteins. This is the
#' internal core that [generate_spm()] wraps: it keeps the per-mutant record of
#' reductions, which [generate_spm()] then reshapes into the lean, model-ready `spm`
#' object. A future motion arm adds its quantities here the same way `dr2_i` / `dr2_n`
#' are added, by calling penm's own per-mutant primitives (`penm::delta_motion_*`,
#' which take `(wt, mut)`) inside this loop, while the mutant is live.
#'
#' The list-column names here keep the index-signature convention (`dr2_ijm`,
#' `dr2_njm`): this is internal, and the letters state the shape exactly, one row per
#' mutant `(j, m)`, each cell a vector over response sites `i` (or modes `n`). The
#' public object renames them to `dr2mat_site` / `dr2mat_mode`; [generate_spm()]
#' is the boundary between the two vocabularies.
#'
#' @param wt Wild-type protein structure with an elastic network model, as
#'   returned by [penm::set_enm()], which also fixes how it mutates.
#' @param pdb_site_active Optional integer vector of active-site residue numbers
#'   (PDB numbering).
#' @return A tibble with one row per mutant `(j, m)`, `m` the letter given to site `j`
#'   (`m = 0` is the wild type),
#'   carrying all measured effects of that mutation: scalar energy changes
#'   (`ddg_dv_jm`, `ddg_tds_jm`, `ddgact_dv_jm`, `ddgact_tds_jm`) and list-columns
#'   `site` / `pdb_site` (the site index to PDB-residue map), `dr2_ijm` (per-site
#'   squared displacement; each cell a `dr2_i` vector over response sites `i` for that
#'   mutant), `mode` (penm's normal-mode index, [penm::get_mode()]), and `dr2_njm`
#'   (per-mode squared contribution; each cell a `dr2_n` vector over response modes `n`).
#' @noRd
generate_spm_core <- function(wt, pdb_site_active = NULL) {
  # Get site information once
  site_vector <- penm::get_site(wt)
  pdb_site_vector <- penm::get_pdb_site(wt)
  # The mode index is penm's own (`nma$mode`), not a locally counted one: asking penm
  # means the label cannot drift from the modes it names if penm's mode set ever changes.
  mode_vector <- penm::get_mode(wt)

  # How the protein mutates was fixed by penm::set_enm(): the mutational model, the
  # ensemble, and the size of the sequence alphabet.
  param <- penm::get_enm_param(wt)
  alphabet_size <- param$alphabet_size

  # Initialize empty results tibble
  results <- tibble()

  # Process each site and letter
  for (j in site_vector) {
    for (m in 0:(alphabet_size - 1)) {
      # Generate single mutant (m=0 is wild type)
      mut <- penm::get_mutant_site(wt, j, m)
      # A genm mutant comes back without normal modes, which the entropy terms below
      # need; an lfenm mutant keeps wt's.
      if (param$mut_model == "genm") mut <- penm::set_enm_nma(mut)

      # Calculate energies
      ddg_dv_jm <- penm::ddg_dv(wt, mut)
      ddg_tds_jm <- penm::ddg_tds(wt, mut)

      # Calculate activity energies if requested
      if (!is.null(pdb_site_active)) {
        ddgact_dv_jm <- penm::ddgact_dv(wt, mut, pdb_site_active = pdb_site_active)
        ddgact_tds_jm <- penm::ddgact_tds(wt, mut, pdb_site_active = pdb_site_active)
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
      if (j %% 50 == 0 && m == alphabet_size - 1) gc()  # Occasional garbage collection
    }
  }

  return(results)
}

#' Run a single-point-mutation scan over a protein
#'
#' Performs a complete single-point mutational scan, giving every residue each other
#' letter of the sequence alphabet and calculating mutational changes of energy and
#' structure.
#'
#' How the protein mutates is set when it is built, by [penm::set_enm()]: the
#' mutational model (`mut_model`), the spread of a letter's shift of a spring length
#' (`mut_dl_sigma`), which springs a mutation leaves alone (`mut_sd_min`), the size of
#' the sequence alphabet (`alphabet_size`) and which sample of mutants is drawn
#' (`ensemble`). The scan reads them from `wt`.
#'
#' Each site carries a letter of the alphabet, numbered `0` to `alphabet_size - 1`; the
#' wild type carries letter `0`, the residue in the pdb, at every site. The scan gives
#' every site `j` each of the other letters in turn, `1` to `alphabet_size - 1`, one
#' mutant per letter. Letters are labels within a site: the same letter at two sites has
#' unrelated effects (see the Details of [penm::set_enm()]). Giving `j` a letter shifts
#' the equilibrium length of every spring touching `j`, each by an independent draw from
#' `Normal(0, mut_dl_sigma)` (only springs between sites at least `mut_sd_min` apart in
#' sequence are shifted). Those altered lengths put the structure under strain, and
#' relaxing it displaces the residues. The scan stores the squared displacements.
#'
#' Experimental, and not covered by the tests: with `mut_model = "genm"`, penm returns
#' each mutant without normal modes; the scan adds them with [penm::set_enm_nma()], since
#' the entropic energy terms need them.
#'
#' The squared displacements are stored in two matrices, `dr2mat_site` and
#' `dr2mat_mode`, both with one row per mutant, in scan order: row `k` of either is the
#' mutant described by row `k` of `energy_data`. They differ in what a column is. In
#' `dr2mat_site` a column is a residue, and cell `[k, i]` is the squared displacement
#' of residue `i` (Å²) under mutant `k`. In `dr2mat_mode` a column is a normal mode,
#' and cell `[k, n]` is the squared projection of that same displacement onto mode `n`.
#' The modes are orthonormal, so each row sums to the same total in both matrices: they
#' hold one displacement written in two bases.
#'
#' @param wt The protein to mutate, carrying its elastic network model and how it
#'   mutates, as returned by [penm::set_enm()]. The scan produces
#'   `n_sites * (alphabet_size - 1)` mutants in all, which is the row count of both
#'   matrices.
#' @param pdb_site_active Active-site residue numbers, in **PDB numbering**. If `NULL`,
#'   no activity energy is computed and `energy_data$ddgact` is `NA` throughout, which
#'   leaves the model's activity term unusable.
#' @return An `spm` object (a classed list) with five elements:
#'   \describe{
#'     \item{`energy_data`}{One row per mutant, aligned with the matrix rows: `j` (the
#'       mutated site), `m` (the letter `j` was given), `ddg` (change in folding
#'       free energy) and `ddgact` (change in activity free energy, or `NA`; see
#'       `pdb_site_active`).}
#'     \item{`dr2mat_site`}{`[mutant x residue]` squared displacements (Å²).}
#'     \item{`dr2mat_mode`}{`[mutant x mode]` squared displacements.}
#'     \item{`site_map`}{The PDB residue number (`pdb_site`) for each `dr2mat_site`
#'       column position (`site`).}
#'     \item{`mode_map`}{The mode index for each `dr2mat_mode` column.}
#'   }
#' @seealso [penm::set_enm()] (builds the elastic network model this needs, and sets
#'   how it mutates); [penm::penm_ensemble] (what `ensemble` means);
#'   [calculate_profiles()] and [fit_lrmsd_msa_site()] (consume the returned object).
#' @family spm
#' @examples
#' if (requireNamespace("bio3d", quietly = TRUE)) {
#'   ex  <- function(f) system.file("extdata", f, package = "msamodel")
#'   wt  <- penm::set_enm(bio3d::read.pdb(ex("1d6o_A.pdb")), node = "ca",
#'                        model = "ming_wall", d_max = 10.5, ensemble = 1L)
#'   act <- read.csv(ex("1d6o_A_active_site.csv"))
#'   spm <- generate_spm(wt, pdb_site_active = act$pdb_site)
#'
#'   dim(spm$dr2mat_site)   # mutants x residues
#'   dim(spm$dr2mat_mode)   # mutants x modes
#'
#'   # The same displacement in two bases: each row sums to the same total.
#'   all.equal(rowSums(spm$dr2mat_site), rowSums(spm$dr2mat_mode))
#' }
#' @export
generate_spm <- function(wt, pdb_site_active = NULL) {
  scan <- generate_spm_core(wt, pdb_site_active = pdb_site_active)

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
